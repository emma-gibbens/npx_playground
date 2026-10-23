%% swr_validation_heatmaps.m
%
% Validation figure: for each trial type (Reward / Airpuff / Neutral),
% produces a two-panel heatmap figure where every row is one neuron:
%
%   Left  : z-scored activity locked to the SPECIFIC STIMULUS the neuron
%            responds to (increase / excitatory direction only):
%              - odor-selective neurons  → locked to odor-1 onset
%              - reward-responsive neurons  → locked to reward delivery
%              - airpuff-responsive neurons → locked to airpuff delivery
%            If a neuron qualifies via both odor and outcome paths for the
%            same trial type, the outcome lock takes precedence (more
%            specific).  Falls back to odor-1 onset if the outcome event
%            timestamps are unavailable in ExpKeys / all_events.mat.
%
%   Right : z-scored activity aligned to ALL SWR onsets in the session
%            (not just participated SWRs), to reveal reactivation strength
%            across the full SWR distribution.
%
% NEURON INCLUSION  (mirrors swr_participation_by_outcome.m / coactivity v2)
%   A neuron enters the heatmap for trial type T if ANY of:
%     - odor_selectivity_eff links to a T-outcome odor
%         AND odor_selectivity_dir_eff == 'increase' / 'excitatory'
%     - T == 'Reward'  AND reward_responsive == 'Yes'
%         AND reward_direction == 'increase'
%     - T == 'Airpuff' AND airpuff_responsive == 'Yes'
%         AND airpuff_direction == 'increase'
%   Neurons with region == 'outside' / 'outside_CA1' or
%   skip_low_baseline_fr == 'Yes' are excluded.
%
% ODOR-SELECTIVE UPGRADE
%   If odor_selectivity == 'None' but odor_responsivity lists exactly ONE
%   odor, that neuron is promoted to odor-selective for that odor.
%   Direction is looked up from py_stats_responsivity_statistic.csv
%   (positive statistic → 'increase'; negative → 'decrease').
%   If odor_responsivity lists 2+ odors the neuron is NOT promoted.
%
% OUTPUT
%   allSessions_heatmap_Reward.png
%   allSessions_heatmap_Airpuff.png
%   allSessions_heatmap_Neutral.png
%   Rows sorted by peak activity time in the 0–0.5 s post-stimulus window.
%   Row labels: "Subject-Date | probeID_unitID"
%
% DEPENDENCIES
%   TSDpeth_fast, getMUA, zscore_tsd, ts, iv, restrict, ResizeIV,
%   SelectTS, load_NPX_spike, LoadExpKeys, FindFiles
%
% EG – June 2026

%% =========================================================================
%% Parameters
%% =========================================================================

folder_list = {'D:\vCA1\M656\preprocessed\M656_2026_03_25'
               'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
               'D:\vCA1\M649\preprocessed\M649_2026_04_01'
               'D:\vCA1\M593\preprocessed\M593_2025_12_26'
               'D:\vCA1\M593\preprocessed\M593_2025_12_28'
               'D:\vCA1\M653\preprocessed\M653_2026_02_04'
               };

% PETH config – stimulus-locked and SWR-locked use identical windows so the
% two heatmap panels are directly comparable. The stimulus panel is a
% validation check: does this neuron actually fire to its tagged stimulus?
% A wide trial-structure window is not needed here.
cfg_peth_trial        = [];
cfg_peth_trial.window = [-1 1];
cfg_peth_trial.dt     = 0.01;
cfg_peth_trial.mode   = 'interp';

% PETH config – SWR-locked (aligned to SWR onset, ALL SWRs)
cfg_peth_swr        = [];
cfg_peth_swr.window = [-0.5 0.5];
cfg_peth_swr.dt     = 0.005;   % 5 ms — resolves SWR-scale dynamics
cfg_peth_swr.mode   = 'interp';

% SWR interval expansion (must match the SWR detection script)
cfg_iv_swr   = [];
cfg_iv_swr.d = [0, 0];

% Colourmap for heatmaps – diverging blue->white->red
n_cmap = 256;
half   = n_cmap / 2;
CMAP   = [linspace(0,1,half)',  linspace(0,1,half)',  ones(half,1); ...
          ones(half,1),         linspace(1,0,half)',  linspace(1,0,half)'];

% Minimum SWRs required to include a session's SWR column
SWR_MIN_COUNT = 20;

%% =========================================================================
%% Accumulators – collect across all sessions before plotting
%% =========================================================================

acc = struct();
for tt = {'Reward','Airpuff','Neutral'}
    tti = tt{1};
    acc.(tti).row_labels = {};   % "Session | label" strings
    acc.(tti).stim_mat   = [];   % [nNeurons × nTimeBins_stim]
    acc.(tti).swr_mat    = [];   % [nNeurons × nTimeBins_swr]
    acc.(tti).tvec_stim  = [];
    acc.(tti).tvec_swr   = [];
    acc.(tti).lock_keys  = {};   % per-row: odor name or 'reward'/'airpuff'
    acc.(tti).lock_types = {};   % per-row: 'odor' or 'outcome'
    acc.(tti).regions    = {};   % per-row: 'dorsal' or 'ventral'
end

%% =========================================================================
%% Main loop – sessions
%% =========================================================================

for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});
    rng(2024);
    LoadExpKeys;

    sessionID = strcat(ExpKeys.subject, '-', ExpKeys.date);
    fprintf('\n=== %s ===\n', sessionID);

    % ------------------------------------------------------------------
    % Load neuron_summary.csv
    % ------------------------------------------------------------------
    csv_path = 'neuron_summary.csv';
    if ~isfile(csv_path)
        warning('neuron_summary.csv not found in %s. Skipping.', folder_list{iFolder});
        continue;
    end
    data = readtable(csv_path, 'VariableNamingRule', 'preserve');

    % ------------------------------------------------------------------
    % Load py_stats_responsivity_statistic.csv (for direction lookup)
    % ------------------------------------------------------------------
    stats_path = 'py_stats_responsivity_statistic.csv';
    if ~isfile(stats_path)
        warning('Stats file not found in %s. Single-odor direction lookup unavailable.', ...
            folder_list{iFolder});
        stats_data = [];
    else
        stats_data = readtable(stats_path, 'VariableNamingRule', 'preserve');
    end

    % ------------------------------------------------------------------
    % Build odor → outcome map from ExpKeys
    % ------------------------------------------------------------------
    odor_to_outcome = containers.Map;
    for i = 1:length(ExpKeys.reward_odors)
        odor_to_outcome(ExpKeys.reward_odors{i}) = 'Reward';
    end
    for i = 1:length(ExpKeys.airpuff_odors)
        odor_to_outcome(ExpKeys.airpuff_odors{i}) = 'Airpuff';
    end
    for i = 1:length(ExpKeys.neutral_odors)
        odor_to_outcome(ExpKeys.neutral_odors{i}) = 'Neutral';
    end

    % ------------------------------------------------------------------
    % Parse neuron_summary rows – mirrors the validated selectivity logic
    % ------------------------------------------------------------------
    % qual_cells: struct array with fields:
    %   .global_id_norm  – label as it appears in S.label (imec0_56 format)
    %   .global_id_raw   – original CSV value (imec0.shank0.56)
    %   .qual_for        – cell array of trial types {'Reward',...}
    %   .lock_to         – cell array (same length as qual_for): for odor-
    %                      locked neurons this is the ACTUAL ODOR NAME (e.g.
    %                      'CS+2') exactly as it appears in evt.label and
    %                      ExpKeys; for outcome-locked neurons it is 'reward'
    %                      or 'airpuff'. This means odor-2-selective neurons
    %                      are locked to odor-2 onset, not odor-1 onset.
    %   .lock_type       – cell array (same length as qual_for): 'odor' or
    %                      'outcome'. Used by the plotting section to choose
    %                      reference lines without relying on odor name strings.
    %   .region  – 'dorsal' or 'ventral' (from neuron_summary)
    qual_cells = struct('global_id_norm',{},'global_id_raw',{},'qual_for',{}, ...
                    'lock_to',{},'lock_type',{},'region',{});

    for row = 1:height(data)

        % Read fields
        global_id_raw   = data.global_id{row};
        global_id_split = split(global_id_raw, '.');
        global_id       = strcat(global_id_split{1}, '_', global_id_split{3});
        % → 'imec0_56' (matches S.label format)

        odor_sel        = data.odor_selectivity{row};
        odor_sel_dir    = data.odor_selectivity_direction{row};
        odor_resp_raw   = data.odor_responsivity{row};
        airpuff_resp    = data.airpuff_responsive{row};
        airpuff_dir     = data.airpuff_direction{row};
        reward_resp     = data.reward_responsive{row};
        reward_dir      = data.reward_direction{row};
        region          = data.region{row};
        skip_fr         = data.skip_low_baseline_fr{row};

        % Exclusion criteria
        if strcmpi(region, 'outside') || strcmpi(region, 'outside_CA1')
            continue;
        end
        if strcmpi(skip_fr, 'Yes')
            continue;
        end

        % Detect empty / NaN / 'None' for odor fields
        odor_sel_is_empty  = isEmptyOrNone(odor_sel);
        odor_resp_is_empty = isEmptyOrNone(odor_resp_raw);

        % ---- Odor-selective upgrade ----------------------------------------
        % Promote single-odor-responsive neurons not already odor-selective.
        if odor_sel_is_empty && ~odor_resp_is_empty
            responsive_odors = strtrim(strsplit(odor_resp_raw, ','));
            responsive_odors = responsive_odors(~cellfun('isempty', responsive_odors));

            if isscalar(responsive_odors)
                odor_sel = responsive_odors{1};

                % Look up response direction from stats CSV
                odor_sel_dir = 'None';
                if ~isempty(stats_data)
                    stats_row_idx = find(strcmp(stats_data.global_unit_id, global_id_raw));
                    if ~isempty(stats_row_idx)
                        stats_cols = stats_data.Properties.VariableNames;
                        col_match  = find(strncmpi(stats_cols, odor_sel, length(odor_sel)));
                        if ~isempty(col_match)
                            stat_val = stats_data{stats_row_idx(1), col_match(1)};
                            if stat_val > 0
                                odor_sel_dir = 'increase';
                            else
                                odor_sel_dir = 'decrease';
                            end
                        else
                            warning('No stats column matching odor "%s" for cell %s', ...
                                odor_sel, global_id_raw);
                        end
                    else
                        warning('Cell %s not found in stats CSV; direction set to None', global_id_raw);
                    end
                end

                fprintf('  Promoted single-odor-responsive cell %s (odor: %s, dir: %s)\n', ...
                    global_id, odor_sel, odor_sel_dir);
            end
            % 2+ responsive odors and not already selective → leave as-is
        end

        % Normalise still-empty odor fields to 'None'
        if isEmptyOrNone(odor_sel)
            odor_sel     = 'None';
            odor_sel_dir = 'None';
        end

        % ---- Classification flags ------------------------------------------
        is_odor_sel = ~strcmpi(odor_sel, 'None') && ~isempty(odor_sel);
        is_reward   = strcmpi(reward_resp,  'Yes');
        is_airpuff  = strcmpi(airpuff_resp, 'Yes');

        if ~is_odor_sel && ~is_reward && ~is_airpuff
            continue;
        end

        % ---- Determine which trial types qualify (increase direction only) --
        qual_for  = {};
        lock_to   = {};   % parallel: odor name (e.g. 'CS+1') or 'reward'/'airpuff'
        lock_type = {};   % parallel: 'odor' or 'outcome'

        % Odor path: only excitatory/increase direction.
        % Store the actual odor identity so the PETH locks to that specific
        % odor's onset, regardless of whether it is odor-1 or odor-2 in the
        % trial structure.
        if is_odor_sel && (strcmpi(odor_sel_dir,'increase') || strcmpi(odor_sel_dir,'excitatory'))
            if isKey(odor_to_outcome, odor_sel)
                this_tt         = odor_to_outcome(odor_sel);
                qual_for{end+1} = this_tt;
                lock_to{end+1}  = odor_sel;    % actual odor name
                lock_type{end+1}= 'odor';
            end
        end

        % Direct outcome paths
        if is_reward  && strcmpi(reward_dir,  'increase')
            qual_for{end+1}  = 'Reward';
            lock_to{end+1}   = 'reward';
            lock_type{end+1} = 'outcome';
        end
        if is_airpuff && strcmpi(airpuff_dir, 'increase')
            qual_for{end+1}  = 'Airpuff';
            lock_to{end+1}   = 'airpuff';
            lock_type{end+1} = 'outcome';
        end

        % De-duplicate: if a neuron qualifies for the same trial type via
        % both odor AND outcome paths, prefer the outcome lock (more specific).
        [qual_for_u, ia] = unique(qual_for, 'last');   % 'last' keeps outcome entry
        qual_for  = qual_for_u;
        lock_to   = lock_to(ia);
        lock_type = lock_type(ia);

        if isempty(qual_for)
            continue;
        end

        n = length(qual_cells) + 1;
        qual_cells(n).global_id_norm = global_id;
        qual_cells(n).global_id_raw  = global_id_raw;
        qual_cells(n).qual_for       = qual_for;
        qual_cells(n).lock_to        = lock_to;
        qual_cells(n).lock_type      = lock_type;
        qual_cells(n).region = region;
    end

    if isempty(qual_cells)
        fprintf('  [skip] No qualifying neurons in %s\n', sessionID);
        clearvars -except folder_list iFolder acc cfg_peth_trial cfg_peth_swr ...
            cfg_iv_swr CMAP SWR_MIN_COUNT;
        continue;
    end

    fprintf('  %d qualifying neurons found.\n', length(qual_cells));

    % ------------------------------------------------------------------
    % Load spikes (imec0 + optional imec1)
    % ------------------------------------------------------------------
    cfg_md          = [];
    record_depth    = ExpKeys.probe1_Depth;
    cfg_md.maxdepth = 1000 * record_depth;

    S = load_NPX_spike('clean_units_imec0.mat');
    S.usr.recording_depths = (cfg_md.maxdepth - S.usr.recording_depths) ...
        .* cos(deg2rad(ExpKeys.probe1_insertion_roll));
    S.usr.probeID = repmat({'imec0'}, 1, length(S.t));

    if isfile('clean_units_imec1.mat')
        record_depth2 = ExpKeys.probe2_Depth;
        S2 = load_NPX_spike('clean_units_imec1.mat');
        S2.usr.recording_depths = (1000 * record_depth2 - S2.usr.recording_depths) ...
            .* cos(deg2rad(ExpKeys.probe2_insertion_roll));
        S2.usr.probeID = repmat({'imec1'}, 1, length(S2.t));

        S.label                = cat(2, S.label,                S2.label);
        S.t                    = cat(2, S.t,                    S2.t);
        S.usr.recording_depths = cat(2, S.usr.recording_depths, S2.usr.recording_depths);
        S.usr.shank_ids        = cat(2, S.usr.shank_ids,        S2.usr.shank_ids);
        S.usr.mean_waveforms   = cat(2, S.usr.mean_waveforms,   S2.usr.mean_waveforms);
        S.usr.probeID          = cat(2, S.usr.probeID,          S2.usr.probeID);
        if isfield(S.usr,'channel_ids') && isfield(S2.usr,'channel_ids')
            S.usr.channel_ids  = cat(2, S.usr.channel_ids, S2.usr.channel_ids);
        end
        clear S2;
    end

    % ------------------------------------------------------------------
    % Load LFP and build MUA time base
    % ------------------------------------------------------------------
    csc      = tsd;
    load('imec0_clean_lfp.mat');
    csc.tvec = imec0.lfp_tvec';
    csc.data = double(imec0.lfp_traces');
    csc.label = cellstr(string(imec0.channel_ids).strip().extractAfter('ap#')');
    csc.cfg.depths = imec0.depths;

    Fs = struct;
    Fs.SamplingFrequency = double(imec0.lfp_fs);
    csc.cfg.hdr{1,1} = Fs;
    csc.cfg.depths = (cfg_md.maxdepth - csc.cfg.depths) ...
        .* cos(deg2rad(ExpKeys.probe1_insertion_roll));

    [~, sort_idx]  = sort(csc.cfg.depths);
    csc.data       = csc.data(sort_idx, :);
    csc.label      = csc.label(sort_idx);
    csc.cfg.depths = csc.cfg.depths(sort_idx);
    clear sort_idx imec0;

    cfg_MUA      = [];
    cfg_MUA.tvec = csc.tvec;
    MUA          = getMUA(cfg_MUA, S);

    % ------------------------------------------------------------------
    % ------------------------------------------------------------------
    % Load events and build event objects
    % ------------------------------------------------------------------
    load('all_events.mat');

    % odor_evts_map: containers.Map  odorName -> ts object (edge-trimmed)
    % Built for every odor in the session so any odor-selective neuron,
    % regardless of whether it prefers odor-1 or odor-2 in the trial
    % structure, is locked to the correct stimulus onset.
    all_session_odors = [ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors];
    odor_evts_map     = containers.Map;
    for iOd = 1:length(all_session_odors)
        od  = all_session_odors{iOd};
        idx = find(strcmp(evt.label, od));
        if isempty(idx)
            fprintf('  [warn] odor "%s" not found in evt.label - skipping.\n', od);
            continue;
        end
        times = getEvtTimes(evt, idx);
        t_evt = ts({times});
        if ~isempty(t_evt.t{1})
            rm = (t_evt.t{1} + cfg_peth_trial.window(1) < MUA.tvec(1)) | ...
                 (t_evt.t{1} + cfg_peth_trial.window(2) > MUA.tvec(end));
            t_evt.t{1}(rm) = [];
        end
        odor_evts_map(od) = t_evt;
    end

    % trial_evts: outcome-locked events (reward / airpuff delivery)
    % Labels are fixed strings in all_events.mat: 'Reward ON' and 'Airpuff ON'
    trial_evts         = struct();
    trial_evts.reward  = [];
    trial_evts.airpuff = [];

    idx_r = find(strcmp(evt.label, 'Reward ON'));
    if ~isempty(idx_r)
        r_times = getEvtTimes(evt, idx_r);
        t_r     = ts({r_times});
        if ~isempty(t_r.t{1})
            rm = (t_r.t{1} + cfg_peth_trial.window(1) < MUA.tvec(1)) | ...
                 (t_r.t{1} + cfg_peth_trial.window(2) > MUA.tvec(end));
            t_r.t{1}(rm) = [];
        end
        trial_evts.reward = t_r;
    else
        fprintf('  [warn] "Reward ON" not found in evt.label for %s\n', sessionID);
    end

    idx_a = find(strcmp(evt.label, 'Airpuff ON'));
    if ~isempty(idx_a)
        a_times = getEvtTimes(evt, idx_a);
        t_a     = ts({a_times});
        if ~isempty(t_a.t{1})
            rm = (t_a.t{1} + cfg_peth_trial.window(1) < MUA.tvec(1)) | ...
                 (t_a.t{1} + cfg_peth_trial.window(2) > MUA.tvec(end));
            t_a.t{1}(rm) = [];
        end
        trial_evts.airpuff = t_a;
    else
        fprintf('  [warn] "Airpuff ON" not found in evt.label for %s\n', sessionID);
    end


    % ------------------------------------------------------------------
    % Load all SWR files; pool onset times per probe
    % ------------------------------------------------------------------
    AllSWR   = FindFiles('*SWR.mat');
    swr_data = struct('probeID',{}, 'location',{}, 'iv',{});

    for currFile = 1:length(AllSWR)
        swrFile = AllSWR{currFile};

        if     contains(swrFile,'imec0'), thisProbe = 'imec0';
        elseif contains(swrFile,'imec1'), thisProbe = 'imec1';
        else,  continue;
        end

        if     contains(swrFile,'dorsal'),  thisLoc = 'dorsal';
        elseif contains(swrFile,'ventral'), thisLoc = 'ventral';
        else
            fprintf('  [skip SWR] cannot determine location: %s\n', swrFile);
            continue;
        end

        loaded = load(swrFile);
        if contains(swrFile,'dorsal')
            rawSWR = loaded.dorsal_swr;
        else
            rawSWR = loaded.ventral_swr;
        end

        if length(rawSWR.iv.tstart) < SWR_MIN_COUNT
            fprintf('  [skip SWR] too few SWRs (%d): %s\n', ...
                length(rawSWR.iv.tstart), swrFile);
            continue;
        end

        n = length(swr_data) + 1;
        swr_data(n).probeID  = thisProbe;
        swr_data(n).location = thisLoc;
        swr_data(n).iv       = ResizeIV(cfg_iv_swr, rawSWR.iv);
    end

    if isempty(swr_data)
        fprintf('  [skip folder] no valid SWR files.\n');
        clearvars -except folder_list iFolder acc cfg_peth_trial cfg_peth_swr ...
            cfg_iv_swr CMAP SWR_MIN_COUNT;
        continue;
    end

    % Pool SWR onsets per probe AND location
    all_swr_onsets_by_probe_loc = struct();
    for thisPr = {'imec0','imec1'}
        pr = thisPr{1};
        for thisLoc = {'dorsal','ventral'}
            loc  = thisLoc{1};
            field = [pr '_' loc];        % 'imec0_dorsal', 'imec0_ventral', etc.
            onsets = [];
            idx_match = find(strcmp({swr_data.probeID}, pr) & ...
                strcmp({swr_data.location}, loc));
            for k = idx_match
                onsets = [onsets; swr_data(k).iv.tstart(:)]; %#ok<AGROW>
            end
            all_swr_onsets_by_probe_loc.(field) = sort(onsets);
        end
    end
    % Build a fast lookup: normalised label → index in S
    S_label_map = containers.Map(S.label, num2cell(1:length(S.label)));

    % ------------------------------------------------------------------
    % Per-qualifying-neuron loop
    % ------------------------------------------------------------------
    for iQ = 1:length(qual_cells)

        qc         = qual_cells(iQ);
        norm_label = qc.global_id_norm;   % 'imec0_56'

        if ~isKey(S_label_map, norm_label)
            fprintf('  [skip] %s not found in S.label\n', norm_label);
            continue;
        end

        iC         = S_label_map(norm_label);
        this_probe = S.usr.probeID{iC};
        this_region   = lower(qc.region);   % 'dorsal' or 'ventral'
        lookup_field  = [this_probe '_' this_region];   % e.g. 'imec0_ventral'
        this_cell  = SelectTS([], S, iC);

        % --- existing: for trial-locked (stimulus on scale of seconds)
        cfg_MUA_trial        = [];
        cfg_MUA_trial.tvec   = csc.tvec;
        cfg_MUA_trial.sigma  = 0.05;   % 50 ms — appropriate for trial-scale
        MUA_cell_trial  = getMUA(cfg_MUA_trial,  this_cell);
        MUAz_trial      = zscore_tsd(MUA_cell_trial);

        % --- new: for SWR-locked (events ~35–150 ms)
        cfg_MUA_swr        = [];
        cfg_MUA_swr.tvec   = csc.tvec;
        cfg_MUA_swr.sigma  = 0.02;   % 20 ms — resolves within-SWR structure
        MUA_cell_swr  = getMUA(cfg_MUA_swr,  this_cell);
        MUAz_swr      = zscore_tsd(MUA_cell_swr);

        % SWR onsets for this cell's probe, edge-trimmed
        all_swr_onsets = all_swr_onsets_by_probe_loc.(lookup_field);
        if ~isempty(all_swr_onsets)
            edge_ok        = (all_swr_onsets + cfg_peth_swr.window(1) >= MUA.tvec(1)) & ...
                             (all_swr_onsets + cfg_peth_swr.window(2) <= MUA.tvec(end));
            all_swr_onsets = all_swr_onsets(edge_ok);
        end

        % SWR PETH (computed once; same row appears in all qualifying heatmaps)
        if ~isempty(all_swr_onsets)
            swr_evt      = ts({all_swr_onsets(:)});
            out_swr      = TSDpeth_fast(cfg_peth_swr, MUAz_swr, swr_evt);
            swr_row      = out_swr.data(:)';
            tvec_swr     = out_swr.tvec;
        else
            n_bins       = round(diff(cfg_peth_swr.window) / cfg_peth_swr.dt) + 1;
            swr_row      = nan(1, n_bins);
            tvec_swr     = linspace(cfg_peth_swr.window(1), cfg_peth_swr.window(2), n_bins);
        end

        row_label = sprintf('%s | %s', sessionID, norm_label);

        % Add to each qualifying trial-type accumulator
        for iT = 1:length(qc.qual_for)
            tti       = qc.qual_for{iT};
            lock_key  = qc.lock_to{iT};    % odor name, 'reward', or 'airpuff'
            ltype     = qc.lock_type{iT};  % 'odor' or 'outcome'

            % Select event timestamps based on lock type
            if strcmp(ltype, 'odor')
                % lock_key is the actual odor name; look it up in odor_evts_map
                if isKey(odor_evts_map, lock_key)
                    lock_evt = odor_evts_map(lock_key);
                else
                    fprintf('  [warn] %s: odor "%s" not in odor_evts_map; skipping.\n', ...
                        norm_label, lock_key);
                    continue;
                end
            else
                % lock_key is 'reward' or 'airpuff'
                lock_evt = trial_evts.(lock_key);
                if isempty(lock_evt) || isempty(lock_evt.t{1})
                    fprintf('  [warn] %s: %s lock unavailable for %s; skipping.\n', ...
                        norm_label, lock_key, tti);
                    continue;
                end
            end

            trial_evt_r = restrict(lock_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
            if isempty(trial_evt_r.t{1})
                continue;
            end

            out_stim    = TSDpeth_fast(cfg_peth_trial, MUAz_trial, trial_evt_r);
            stim_row    = out_stim.data(:)';
            tvec_stim   = out_stim.tvec;

            acc.(tti).row_labels{end+1} = row_label;
            acc.(tti).stim_mat          = [acc.(tti).stim_mat;  stim_row];
            acc.(tti).swr_mat           = [acc.(tti).swr_mat;   swr_row];
            acc.(tti).tvec_stim         = tvec_stim;
            acc.(tti).tvec_swr          = tvec_swr;
            acc.(tti).lock_keys{end+1}  = lock_key;   % odor name or 'reward'/'airpuff'
            acc.(tti).lock_types{end+1} = ltype;      % 'odor' or 'outcome'
            acc.(tti).regions{end+1}    = lower(qc.region);  % 'dorsal' or 'ventral'
        end

    end % qualifying-neuron loop

    clearvars -except folder_list iFolder acc cfg_peth_trial cfg_peth_swr ...
        cfg_iv_swr CMAP SWR_MIN_COUNT;

end % folder loop

%% =========================================================================
%% Plotting – one figure per trial type
%% =========================================================================

trial_types  = {'Reward',  'Airpuff', 'Neutral'};

for iTT = 1:length(trial_types)
    tti  = trial_types{iTT};
    dat  = acc.(tti);
    nNeur = size(dat.stim_mat, 1);

    if nNeur == 0
        fprintf('[skip] No neurons qualify for %s heatmap.\n', tti);
        continue;
    end

    fprintf('\nPlotting %s heatmap: %d neurons.\n', tti, nNeur);

    % Sort rows: (1) region — dorsal first, ventral second;
    %            (2) within each region, by peak activity time in [0, +1.0 s]
    stim_tvec   = dat.tvec_stim;
    peak_window = stim_tvec >= 0 & stim_tvec <= 1.0;

    % Assign numeric region codes: dorsal = 1, ventral = 2, other = 3
    region_code = ones(nNeur, 1) * 3;
    for iR = 1:nNeur
        if strcmpi(dat.regions{iR}, 'dorsal'),  region_code(iR) = 1;
        elseif strcmpi(dat.regions{iR}, 'ventral'), region_code(iR) = 2;
        end
    end

    % Peak bin within window (NaN rows → pushed to end within group)
    if sum(peak_window) > 0
        sub_mat    = dat.stim_mat(:, peak_window);
        peak_bins  = nan(nNeur, 1);
        for iR = 1:nNeur
            row_vals = sub_mat(iR, :);
            if any(~isnan(row_vals))
                [~, pk]       = max(row_vals);
                peak_bins(iR) = pk;
            else
                peak_bins(iR) = Inf;
            end
        end
    else
        peak_bins = (1:nNeur)';
    end

    % Two-key sort: region code first, peak bin second
    [~, sort_order] = sortrows([region_code, peak_bins]);

    stim_sorted       = dat.stim_mat(sort_order, :);
    swr_sorted        = dat.swr_mat(sort_order,  :);
    labels_sorted     = dat.row_labels(sort_order);
    locks_sorted      = dat.lock_keys(sort_order);   % odor name or 'reward'/'airpuff'
    locktypes_sorted  = dat.lock_types(sort_order);  % 'odor' or 'outcome'
    regions_sorted    = dat.regions(sort_order);

    % Row index of last dorsal neuron (for divider line between regions)
    last_dorsal = sum(region_code(sort_order) == 1);

    % Robust colour limits (1st-99th percentile, both panels share scale)
    % Separate robust limits per panel
    stim_vals = stim_sorted(~isnan(stim_sorted(:)));
    swr_vals  = swr_sorted(~isnan(swr_sorted(:)));

    clim_stim_lo = prctile(stim_vals, 1);
    clim_stim_hi = prctile(stim_vals, 97);
    clim_swr_lo  = prctile(swr_vals,  1);
    clim_swr_hi  = prctile(swr_vals,  97);

    tick_fs = max(5, min(10, floor(200 / nNeur)));

    fig_h = figure('Units','normalized','Position',[0.02 0.05 0.90 0.90]);

    % --- LEFT: stimulus-locked ---
    ax1 = subplot(1, 2, 1);
    imagesc(ax1, stim_tvec, 1:nNeur, stim_sorted);
    clim(ax1, [clim_stim_lo clim_stim_hi]);
    colormap(ax1, CMAP);
    hold(ax1, 'on');
    xline(ax1, 0, '--', 'Color',[0 0 0], 'LineWidth',1.5, 'Label','Stimulus onset');

    % Region divider and labels
    n_dorsal  = last_dorsal;
    n_ventral = nNeur - last_dorsal;
    if n_dorsal > 0 && n_ventral > 0
        yline(ax1, last_dorsal + 0.5, '-', 'Color',[0.2 0.2 0.2], 'LineWidth',2);
    end
    x_lbl = stim_tvec(1) + 0.02 * diff([stim_tvec(1) stim_tvec(end)]);
    if n_dorsal > 0
        text(ax1, x_lbl, n_dorsal/2, ...
            sprintf('Dorsal (n=%d)', n_dorsal), ...
            'Color','w','FontSize',9,'FontWeight','bold', ...
            'VerticalAlignment','middle','HorizontalAlignment','left', ...
            'BackgroundColor',[0 0 0 0.3],'Margin',1);
    end
    if n_ventral > 0
        text(ax1, x_lbl, last_dorsal + n_ventral/2, ...
            sprintf('Ventral (n=%d)', n_ventral), ...
            'Color','w','FontSize',9,'FontWeight','bold', ...
            'VerticalAlignment','middle','HorizontalAlignment','left', ...
            'BackgroundColor',[0 0 0 0.3],'Margin',1);
    end

    % Build x-axis label from the lock types actually present
    has_odor_lock    = any(strcmp(locktypes_sorted, 'odor'));
    has_reward_lock  = any(strcmp(locks_sorted,     'reward'));
    has_airpuff_lock = any(strcmp(locks_sorted,     'airpuff'));
    lock_descs = {};
    if has_odor_lock,    lock_descs{end+1} = 'odor onset'; end
    if has_reward_lock,  lock_descs{end+1} = 'reward delivery'; end
    if has_airpuff_lock, lock_descs{end+1} = 'airpuff delivery'; end
    xlab_str = sprintf('Time from %s (s)', strjoin(lock_descs, ' / '));

    yticks(ax1, 1:nNeur);
    yticklabels(ax1, labels_sorted);
    ax1.TickLabelInterpreter = 'none';
    ax1.FontSize = tick_fs;
    xlabel(ax1, xlab_str, 'FontSize', 12);
    ylabel(ax1, 'Neuron  (session | probe_unit)', 'FontSize', 11);
    title(ax1, sprintf('%s - Stimulus-locked (n=%d)', tti, nNeur), 'FontSize', 13);
    cb1 = colorbar(ax1); cb1.Label.String = 'Z-score';
    axis(ax1, 'tight');

    % --- RIGHT: SWR-locked ---
    ax2 = subplot(1, 2, 2);
    imagesc(ax2, dat.tvec_swr, 1:nNeur, swr_sorted);
    clim(ax2, [clim_swr_lo  clim_swr_hi]);
    colormap(ax2, CMAP);
    hold(ax2, 'on');

    xline(ax2, 0, '--', 'Color',[0 0 0], 'LineWidth',1.5, 'Label','SWR onset');

    % Region divider (same boundary as left panel)
    if n_dorsal > 0 && n_ventral > 0
        yline(ax2, last_dorsal + 0.5, '-', 'Color',[0.2 0.2 0.2], 'LineWidth',2);
    end
    x_lbl2 = dat.tvec_swr(1) + 0.02 * diff([dat.tvec_swr(1) dat.tvec_swr(end)]);
    if n_dorsal > 0
        text(ax2, x_lbl2, n_dorsal/2, ...
            sprintf('Dorsal (n=%d)', n_dorsal), ...
            'Color','w','FontSize',9,'FontWeight','bold', ...
            'VerticalAlignment','middle','HorizontalAlignment','left', ...
            'BackgroundColor',[0 0 0 0.3],'Margin',1);
    end
    if n_ventral > 0
        text(ax2, x_lbl2, last_dorsal + n_ventral/2, ...
            sprintf('Ventral (n=%d)', n_ventral), ...
            'Color','w','FontSize',9,'FontWeight','bold', ...
            'VerticalAlignment','middle','HorizontalAlignment','left', ...
            'BackgroundColor',[0 0 0 0.3],'Margin',1);
    end

    yticks(ax2, 1:nNeur);
    yticklabels(ax2, labels_sorted);
    ax2.TickLabelInterpreter = 'none';
    ax2.FontSize = tick_fs;
    xlabel(ax2, 'Time from SWR onset (s)', 'FontSize', 12);
    ylabel(ax2, 'Neuron  (session | probe_unit)', 'FontSize', 11);
    title(ax2, sprintf('%s – SWR-locked, all SWRs (n=%d)', tti, nNeur), 'FontSize', 13);
    cb2 = colorbar(ax2); cb2.Label.String = 'Z-score';
    axis(ax2, 'tight');

    linkaxes([ax1, ax2], 'y');
    sgtitle(sprintf('Validation heatmaps – %s neurons  (dorsal: %d | ventral: %d)', ...
        tti, n_dorsal, n_ventral), 'FontSize', 14);

    out_fname = sprintf('allSessions_heatmap_%s', tti);
    print(fig_h, '-dpng', '-r150', out_fname);
    fprintf('  Saved: %s.png\n', out_fname);
end

fprintf('\nDone.\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function result = isEmptyOrNone(val)
% Returns true if val is empty, numeric NaN, or the string 'None'.
    result = isempty(val) || strcmpi(val, 'None') || ...
             (isscalar(val) && isnumeric(val) && isnan(val));
end

function times = getEvtTimes(evt, idx)
% Return event times at position idx in evt, or [] if absent.
    if ~isempty(idx)
        times = evt.t{idx};
    else
        times = [];
    end
end