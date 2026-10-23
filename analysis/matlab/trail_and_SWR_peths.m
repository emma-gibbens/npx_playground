%% trial_and_SWR_peths.m
% Generates per-neuron figures containing:
%   (1) A full-trial PETH aligned to odor 1 onset, grouped by outcome category
%       (reward = blue, airpuff = red, neutral = black)
%   (2) A SWR PETH aligned to the onset of SWRs in which the cell fired at
%       least one spike ("participated" SWRs only)
%
% OVERVIEW
%   For each session folder the script:
%     1. Loads spikes (imec0, optionally imec1), LFP, events, and ExpKeys
%        exactly as in peths_task3_odors_catch_SWR.m
%     2. Finds every SWR file in the folder and, for each probe, identifies
%        the SWRs the cell actually fired in (spike_count > 0 inside the
%        expanded SWR window), following the participation logic in
%        odor_selective_swr_participation_dvHC_v4.m
%     3. For each cell plots:
%         - Left panel  : z-scored full-trial PETH (reward / airpuff / neutral)
%         - Right panel : z-scored SWR PETH (participated SWRs only, per probe)
%     4. Saves one PNG per cell
%
% FILTERING
%   Only neurons present in neuron_summary.csv that meet at least one of:
%     - odor_responsivity is not "None"  (odor-responsive)
%     - airpuff_responsive  == "Yes"
%     - reward_responsive   == "Yes"
%   are plotted. All others are skipped.
%
% TRIAL TIMING (from odor 1 onset)
%   Odor 1 onset   :  0.00 s    Odor 1 offset :  1.00 s
%   Odor 2 onset   :  1.50 s    Odor 2 offset :  2.50 s
%   Airpuff        :  2.75 s    Reward        :  3.00 s
%
% DEPENDENCIES
%   mvdmlab / odor-pixels codebase (TSDpeth_fast, getMUA, zscore_tsd,
%   SelectTS, ts, iv, restrict, ResizeIV, MakeQfromS, load_NPX_spike,
%   LoadExpKeys, FindFiles, ccf)
%
% INPUT FILES (per session folder)
%   clean_units_imec0.mat   – required
%   clean_units_imec1.mat   – optional (dual-probe)
%   imec0_clean_lfp.mat
%   all_events.mat
%   ExpKeys  (via LoadExpKeys)
%   *SWR*.mat files         – one or more, named with 'imec0'/'imec1' and
%                             'dorsal'/'ventral'
%
% OUTPUT
%   One PNG per cell, named:
%   [subject]_[date]_[unit_label]_trial_SWR_peth.png

%% =========================================================================
%% Parameters – edit here
%% =========================================================================

folder_list = {'D:\vCA1\M656\preprocessed\M656_2026_03_25'};

% Path to the neuron summary CSV containing responsivity labels
neuron_summary_path = 'D:\vCA1\M656\preprocessed\M656_2026_03_25\neuron_summary.csv'; 

% SWR window expansion applied to each SWR interval before participation
% check (must match what was used when the SWR files were created)
cfg_iv_swr     = [];
cfg_iv_swr.d   = [-0.025, 0.025];   % +/- 25 ms on each side

% Dorsal/ventral depth thresholds (microns from brain surface, per probe)
DH_threshold_imec0 = 3000;
DH_threshold_imec1 = 3000;

% Full-trial PETH config (aligned to odor 1 onset)
cfg_peth_trial        = [];
cfg_peth_trial.window = [-1 4];
cfg_peth_trial.dt     = 0.01;
cfg_peth_trial.mode   = 'interp';

% SWR PETH config (aligned to SWR onset)
cfg_peth_swr        = [];
cfg_peth_swr.window = [-0.5 0.5];   % 500 ms either side of SWR onset
cfg_peth_swr.dt     = 0.01;
cfg_peth_swr.mode   = 'interp';

% Fixed US latencies from odor 1 onset (for reference lines)
airpuff_latency_from_odor1 = 2.75;
reward_latency_from_odor1  = 3.0;

% Colors
col_reward  = [0 0 1];   % blue
col_airpuff = [1 0 0];   % red
col_neutral = [0 0 0];   % black
col_swr     = [0.4 0 0.6]; % purple – participated SWRs

%% =========================================================================
%% Load neuron summary and build set of responsive unit labels
%% =========================================================================

neuron_tbl = readtable(neuron_summary_path, 'TextType', 'string');

% A neuron passes the filter if ANY of the three criteria are true:
%   1. odor_responsivity is not "None"
%   2. airpuff_responsive == "Yes"
%   3. reward_responsive  == "Yes"
odor_selective    = neuron_tbl.odor_selectivity ~= "None";
airpuff_responsive = neuron_tbl.airpuff_responsive == "Yes";
reward_responsive  = neuron_tbl.reward_responsive  == "Yes";

passes_filter = odor_selective | airpuff_responsive | reward_responsive;

% Store the global_ids of all neurons that pass as a set for fast lookup
responsive_ids = neuron_tbl.global_id(passes_filter);
N = length(responsive_ids);
responsive_ids_fixed = strings(N, 1);
for i = 1:N
    gid = responsive_ids(i);
    short_id = regexprep(gid, 'imec0\.shank\d+\.', 'imec0_'); % I apologize for how unreadable this is
    responsive_ids_fixed(i) = short_id;
end


fprintf('Neuron filter: %d / %d neurons pass (odor- or outcome-responsive).\n', ...
    sum(passes_filter), height(neuron_tbl));

%% =========================================================================
%% Main loop
%% =========================================================================

for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});
    rng(2024);
    LoadExpKeys;

    % ------------------------------------------------------------------
    % Depth calibration
    % ------------------------------------------------------------------
    cfg_md          = [];
    record_depth    = ExpKeys.probe1_Depth;
    cfg_md.maxdepth = 1000 * record_depth;

    % ------------------------------------------------------------------
    % Load spikes (imec0, + imec1 if present)
    % ------------------------------------------------------------------
    S = load_NPX_spike('clean_units_imec0.mat');
    S.usr.recording_depths = (cfg_md.maxdepth - S.usr.recording_depths) ...
        * cos(deg2rad(ExpKeys.probe1_insertion_roll));

    % Track which probe each unit came from (needed to match SWR files)
    S.usr.probeID = repmat({'imec0'}, 1, length(S.t));

    if isfile('clean_units_imec1.mat')
        record_depth2  = ExpKeys.probe2_Depth;
        S2 = load_NPX_spike('clean_units_imec1.mat');
        S2.usr.recording_depths = (1000 * record_depth2 - S2.usr.recording_depths) ...
            * cos(deg2rad(ExpKeys.probe2_insertion_roll));
        S2_probeIDs = repmat({'imec1'}, 1, length(S2.t));

        S.label  = cat(2, S.label,  S2.label);
        S.t      = cat(2, S.t,      S2.t);
        S.usr.recording_depths = cat(2, S.usr.recording_depths, S2.usr.recording_depths);
        S.usr.shank_ids        = cat(2, S.usr.shank_ids,        S2.usr.shank_ids);
        S.usr.mean_waveforms   = cat(2, S.usr.mean_waveforms,   S2.usr.mean_waveforms);
        S.usr.probeID          = cat(2, S.usr.probeID,          S2_probeIDs);
        if isfield(S.usr,'channel_ids') && isfield(S2.usr,'channel_ids')
            S.usr.channel_ids = cat(2, S.usr.channel_ids, S2.usr.channel_ids);
        end
        clear S2 S2_probeIDs
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
        * cos(deg2rad(ExpKeys.probe1_insertion_roll));

    [~, sort_idx] = sort(csc.cfg.depths);
    csc.data       = csc.data(sort_idx,:);
    csc.label      = csc.label(sort_idx);
    csc.cfg.depths = csc.cfg.depths(sort_idx);
    clear sort_idx

    cfg_MUA      = [];
    cfg_MUA.tvec = csc.tvec;
    MUA  = getMUA(cfg_MUA, S);
    MUAz = zscore_tsd(MUA);

    % ------------------------------------------------------------------
    % Load events and identify odor-1 onset times per outcome category
    % ------------------------------------------------------------------
    load('all_events.mat');

    reward_odor1_label  = ExpKeys.reward_odors{1};
    airpuff_odor1_label = ExpKeys.airpuff_odors{1};
    neutral_odor1_label = ExpKeys.neutral_odors{1};

    % Human-readable odor label strings for the figure title
    reward_odors_labels  = strjoin(extractAfter(ExpKeys.reward_odors,  'Odor '), '');
    airpuff_odors_labels = strjoin(extractAfter(ExpKeys.airpuff_odors, 'Odor '), '');
    neutral_odors_labels = strjoin(extractAfter(ExpKeys.neutral_odors, 'Odor '), '');

    reward_odor1_idx  = find(strcmp(evt.label, reward_odor1_label));
    airpuff_odor1_idx = find(strcmp(evt.label, airpuff_odor1_label));
    neutral_odor1_idx = find(strcmp(evt.label, neutral_odor1_label));

    reward_odor1_times  = getEvtTimes(evt, reward_odor1_idx);
    airpuff_odor1_times = getEvtTimes(evt, airpuff_odor1_idx);
    neutral_odor1_times = getEvtTimes(evt, neutral_odor1_idx);

    % ------------------------------------------------------------------
    % Build trial event objects and strip edges
    % ------------------------------------------------------------------
    reward_trial_evt  = ts({reward_odor1_times});
    airpuff_trial_evt = ts({airpuff_odor1_times});
    neutral_trial_evt = ts({neutral_odor1_times});

    trial_evts = {reward_trial_evt, airpuff_trial_evt, neutral_trial_evt};
    for iE = 1:length(trial_evts)
        if ~isempty(trial_evts{iE}.t{1})
            tstart_chk = trial_evts{iE}.t{1} + cfg_peth_trial.window(1);
            tend_chk   = trial_evts{iE}.t{1} + cfg_peth_trial.window(2);
            rm         = tstart_chk < MUA.tvec(1) | tend_chk > MUA.tvec(end);
            trial_evts{iE}.t{1}(rm) = [];
        end
    end
    reward_trial_evt  = trial_evts{1};
    airpuff_trial_evt = trial_evts{2};
    neutral_trial_evt = trial_evts{3};

    % ------------------------------------------------------------------
    % Load all SWR files and collect per-probe SWR intervals
    % ------------------------------------------------------------------
    % swr_data is a struct array with fields:
    %   .probeID   – 'imec0' or 'imec1'
    %   .location  – 'dorsal' or 'ventral'
    %   .iv        – interval struct with tstart / tend (already resized)
    % We keep all SWR sets found; during the cell loop we match each cell
    % to the SWR sets from its own probe.
    AllSWR   = FindFiles('*SWR.mat');
    swr_data = struct('probeID', {}, 'location', {}, 'iv', {});

    for currFile = 1:length(AllSWR)
        swrFile = AllSWR{currFile};

        if contains(swrFile, 'imec0')
            thisProbe = 'imec0';
        elseif contains(swrFile, 'imec1')
            thisProbe = 'imec1';
        else
            continue;
        end

        if contains(swrFile, 'dorsal')
            thisLoc = 'dorsal';
        elseif contains(swrFile, 'ventral')
            thisLoc = 'ventral';
        else
            fprintf('  [skip SWR file] cannot determine dorsal/ventral: %s\n', swrFile);
            continue;
        end

        loaded = load(swrFile);
        if contains(swrFile, 'dorsal')
            rawSWR = loaded.dorsal_swr;
        else
            rawSWR = loaded.ventral_swr;
        end

        % Expand each SWR interval by cfg_iv_swr.d on each side
        expandedIV = ResizeIV(cfg_iv_swr, rawSWR.iv);

        entry.probeID  = thisProbe;
        entry.location = thisLoc;
        entry.iv       = expandedIV;
        swr_data(end+1) = entry; %#ok<AGROW>
    end

    if isempty(swr_data)
        fprintf('  [skip folder] No valid SWR files found in %s\n', folder_list{iFolder});
        clearvars -except folder_list iFolder responsive_ids responsive_ids_fixed neuron_summary_path
        continue;
    end

    % ------------------------------------------------------------------
    % Per-cell loop
    % ------------------------------------------------------------------
    for iC = 1:length(S.t)

        this_cell  = SelectTS([], S, iC);
        this_probe = S.usr.probeID{iC};
        currCell_t = S.t{iC};

        % --------------------------------------------------------------
        % FILTER: skip this neuron if it is not odor- or outcome-responsive
        % The global_id in the CSV is formatted as e.g. "imec0.shank0.14".
        % S.label entries are expected to match this format directly.
        % --------------------------------------------------------------
        cell_label = this_cell.label{1};   % e.g. 'imec0.shank0.14'
        if ~any(strcmp(responsive_ids_fixed, cell_label))
            fprintf('  [skip] %s – not odor/outcome responsive\n', cell_label);
            continue;
        end

        MUA_cell  = getMUA(cfg_MUA, this_cell);
        MUAz_cell = zscore_tsd(MUA_cell);

        % Find which SWR sets belong to this cell's probe
        probe_swr_mask = strcmp({swr_data.probeID}, this_probe);
        probe_swr_sets = swr_data(probe_swr_mask);

        % Collect onset times of SWRs this cell *participated* in,
        % pooled across all SWR sets from the matching probe.
        % A cell "participated" if it fired >= 1 spike inside the SWR window.
        participated_swr_onsets = [];

        for iSet = 1:length(probe_swr_sets)
            this_iv   = probe_swr_sets(iSet).iv;
            nSWR      = length(this_iv.tstart);

            for iSWR = 1:nSWR
                t_start     = this_iv.tstart(iSWR);
                t_end       = this_iv.tend(iSWR);
                spike_count = sum(currCell_t > t_start & currCell_t < t_end);
                if spike_count > 0
                    participated_swr_onsets(end+1) = t_start; %#ok<AGROW>
                end
            end
        end

        % Build the figure: 1 row, 2 panels
        this_fig = figure('WindowState', 'maximized');

        % ==============================================================
        % LEFT PANEL – z-scored full-trial PETH (aligned to odor 1)
        % ==============================================================
        ax_trial = subplot(2, 1, 1);
        hold(ax_trial, 'on');

        leg_labels_trial = {};
        min_ylim_t =  10;
        max_ylim_t = -10;

        % Reward trials
        rwd_evt_r = restrict(reward_trial_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(rwd_evt_r.t{1})
            out_z = TSDpeth_fast(cfg_peth_trial, MUAz_cell, rwd_evt_r);
            plot(ax_trial, out_z.tvec, out_z.data, 'Color', col_reward, 'LineWidth', 2);
            leg_labels_trial{end+1} = sprintf('Reward (n=%d)', length(rwd_evt_r.t{1}));
            min_ylim_t = min(min_ylim_t, min(out_z.data) - 0.1);
            max_ylim_t = max(max_ylim_t, max(out_z.data) + 0.1);
        end

        % Airpuff trials
        apf_evt_r = restrict(airpuff_trial_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(apf_evt_r.t{1})
            out_z = TSDpeth_fast(cfg_peth_trial, MUAz_cell, apf_evt_r);
            plot(ax_trial, out_z.tvec, out_z.data, 'Color', col_airpuff, 'LineWidth', 2);
            leg_labels_trial{end+1} = sprintf('Airpuff (n=%d)', length(apf_evt_r.t{1}));
            min_ylim_t = min(min_ylim_t, min(out_z.data) - 0.1);
            max_ylim_t = max(max_ylim_t, max(out_z.data) + 0.1);
        end

        % Neutral trials
        neu_evt_r = restrict(neutral_trial_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(neu_evt_r.t{1})
            out_z = TSDpeth_fast(cfg_peth_trial, MUAz_cell, neu_evt_r);
            plot(ax_trial, out_z.tvec, out_z.data, 'Color', col_neutral, 'LineWidth', 2);
            leg_labels_trial{end+1} = sprintf('Neutral (n=%d)', length(neu_evt_r.t{1}));
            min_ylim_t = min(min_ylim_t, min(out_z.data) - 0.1);
            max_ylim_t = max(max_ylim_t, max(out_z.data) + 0.1);
        end

        % Reference lines for trial structure
        xline(ax_trial, 0,    '-',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', 'Odor 1 on');
        xline(ax_trial, 1.0,  '-',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'Label', 'Odor 1 off');
        xline(ax_trial, 1.5,  '-',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', 'Odor 2 on');
        xline(ax_trial, 2.5,  '-',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'Label', 'Odor 2 off');
        xline(ax_trial, airpuff_latency_from_odor1, '--', 'Color', col_airpuff, 'LineWidth', 1.5, 'Label', 'Airpuff');
        xline(ax_trial, reward_latency_from_odor1,  '--', 'Color', col_reward,  'LineWidth', 1.5, 'Label', 'Reward');

        ax_trial.YLim      = [min(min_ylim_t, -0.5), max(max_ylim_t, 1)];
        ax_trial.Title.String   = 'Z-scored Trial PETH (aligned to Odor 1 onset)';
        ax_trial.Title.FontSize = 14;
        ax_trial.YLabel.String  = 'Z-score';
        ax_trial.YLabel.FontSize = 14;
        ax_trial.XLabel.String  = 'Time from Odor 1 onset (s)';
        ax_trial.XLabel.FontSize = 14;
        legend(ax_trial, leg_labels_trial, 'Location', 'best', 'FontSize', 10, 'Box', 'off');

        % ==============================================================
        % RIGHT PANEL – z-scored SWR PETH (participated SWRs only)
        % ==============================================================
        ax_swr = subplot(2, 1, 2);
        hold(ax_swr, 'on');

        if ~isempty(participated_swr_onsets)
            % Remove SWR onsets that are too close to recording edges
            swr_tstart_chk = participated_swr_onsets + cfg_peth_swr.window(1);
            swr_tend_chk   = participated_swr_onsets + cfg_peth_swr.window(2);
            edge_ok        = swr_tstart_chk >= MUA.tvec(1) & swr_tend_chk <= MUA.tvec(end);
            participated_swr_onsets_trimmed = participated_swr_onsets(edge_ok);

            if ~isempty(participated_swr_onsets_trimmed)
                swr_evt = ts({participated_swr_onsets_trimmed(:)});
                out_swr_z = TSDpeth_fast(cfg_peth_swr, MUAz_cell, swr_evt);

                plot(ax_swr, out_swr_z.tvec, out_swr_z.data, ...
                    'Color', col_swr, 'LineWidth', 2);

                swr_leg = sprintf('Participated SWRs (n=%d)', ...
                    length(participated_swr_onsets_trimmed));
            else
                swr_leg = 'No valid participated SWRs';
            end
        else
            swr_leg = 'No participated SWRs';
        end

        % Reference line at SWR onset
        xline(ax_swr, 0, '--', 'Color', [0 0 0], 'LineWidth', 1.5, 'Label', 'SWR onset');

        ax_swr.Title.String    = 'Z-scored SWR PETH (participated SWRs only)';
        ax_swr.Title.FontSize  = 14;
        ax_swr.YLabel.String   = 'Z-score';
        ax_swr.YLabel.FontSize = 14;
        ax_swr.XLabel.String   = 'Time from SWR onset (s)';
        ax_swr.XLabel.FontSize = 14;
        legend(ax_swr, {swr_leg}, 'Location', 'best', 'FontSize', 10, 'Box', 'off');

        % ==============================================================
        % Figure title and save
        % ==============================================================
        fname = strcat(ExpKeys.subject, '_', strrep(ExpKeys.date, '-', '_'), ...
            '_', this_cell.label{1}, '_trial_SWR_peth');
        tname = sprintf('%s  |  Probe: %s  |  Reward: %s  Airpuff: %s  Neutral: %s', ...
            fname, this_probe, reward_odors_labels, airpuff_odors_labels, neutral_odors_labels);
        sgtitle(tname, 'Interpreter', 'none', 'FontSize', 11);

        print(this_fig, '-dpng', fname);
        close(this_fig);

    end % cell loop

    clearvars -except folder_list iFolder responsive_ids responsive_ids_fixed neuron_summary_path

end % folder loop


%% =========================================================================
%% LOCAL HELPER
%% =========================================================================

function times = getEvtTimes(evt, idx)
% Returns spike times for the event at position idx in evt, or [] if absent.
    if ~isempty(idx)
        times = evt.t{idx};
    else
        times = [];
    end
end