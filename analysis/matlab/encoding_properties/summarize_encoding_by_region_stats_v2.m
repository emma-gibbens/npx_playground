%% Summarize odor, reward, and airpuff encoding by region
% Reports: for each region (dorsal / ventral), what % of CA1 neurons
% are odor-selective, reward-responsive, and airpuff-responsive.
%
% Denominator = CA1 neurons loaded from clean_units_imec0/1.mat,
% after depth correction, dorsal/ventral assignment, AND CA1 depth-band
% filtering (each neuron must be within ±CA1_DEPTH_THRESH_UM of the best
% SWR channel on its probe/region). This exactly matches the cell
% population used in swr_participation_stimulus_encoding_v4.m.
%
% The CA1 filter runs per SWR file (one per probe×region combination),
% exactly mirroring v4. A neuron is kept in the denominator if it passes
% the CA1 band filter for at least one SWR file on its probe.
%
% Encoding numerators come from neuron_summary.csv, using the same
% inclusion and single-odor-promotion logic as generateSelectiveTS_v3.m.
% 
% NEW TO VERSION 2: ENCODING OVERLAP BINS 
%   This file tries to account for the fact that odor-selective,
%   reward-responsive, and airpuff-responsive classifications may overlap.
%   So there are 8 mutually exclusive bins (odor only, reward only, airpuff
%   only, odor + reward, odor + airpuff, reward + airpuff, all three, none)
%
% Denominator: ALL surviving (post-decorr-mask, post-CA1-band-filter)
%   clean_units neurons assigned to dorsal/ventral, regardless of whether
%   they appear in neuron_summary.csv. Neurons in the CSV that are flagged
%   outside or skip_low_baseline_fr are excluded from encoding numerators
%   but still count in the total.
%
% Decorrelation mask: neurons on shanks with a non-negative decorr time
%   in ExpKeys are excluded from BOTH the denominator and the encoding
%   numerators (treated as if never recorded). Mirrors v4.
%
% Output: per-session table + aggregate totals across all sessions.
%
% EG June 2026

function summarize_encoding_by_region_stats()

%% -------------------------------------------------------------------------
%% Parameters  (match swr_participation_stimulus_encoding_v4.m)
%% -------------------------------------------------------------------------

folder_list = {
    'D:\vCA1\M656\preprocessed\M656_2026_03_25'
    'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
    'D:\vCA1\M649\preprocessed\M649_2026_04_01'
    'D:\vCA1\M593\preprocessed\M593_2025_12_26'
    'D:\vCA1\M593\preprocessed\M593_2025_12_28'
    'D:\vCA1\M653\preprocessed\M653_2026_02_04'
};

% folder_list = {
%     'D:\vCA1\M560\preprocessed\M560_2025_12_19'
%     'D:\vCA1\M593\preprocessed\M593_2025_12_26'
%     'D:\vCA1\M593\preprocessed\M593_2025_12_28'
%     'D:\vCA1\M653\preprocessed\M653_2026_02_04'
%     'D:\vCA1\M656\preprocessed\M656_2026_03_25'
%     'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
%     'D:\vCA1\M649\preprocessed\M649_2026_04_01'
%     'D:\vCA1\M656\preprocessed\M656_2026_03_26_conservativeSpikeSorting'
%     };

SWR_MIN_COUNT = 20;   % min SWRs required to use a SWR file (mirrors v4)

%% -------------------------------------------------------------------------
%% Accumulators across all sessions
%% -------------------------------------------------------------------------

regions = {'Dorsal', 'Ventral'};
outcomes = {'Reward', 'Airpuff', 'Neutral'};

agg = struct();
for r = 1:length(regions)
    rn = regions{r};
    agg.(rn).total        = 0;
    agg.(rn).odor         = 0;
    agg.(rn).reward       = 0;
    agg.(rn).airpuff      = 0;
    % Odor-selective neurons broken down by outcome category
    agg.(rn).odor_reward  = 0;
    agg.(rn).odor_airpuff = 0;
    agg.(rn).odor_neutral = 0;
    % Overlap counts (neurons positive for >1 encoding type)
    agg.(rn).odor_only        = 0;
    agg.(rn).reward_only      = 0;
    agg.(rn).airpuff_only     = 0;
    agg.(rn).odor_and_reward  = 0;   % odor-selective AND reward-responsive
    agg.(rn).odor_and_airpuff = 0;   % odor-selective AND airpuff-responsive
    agg.(rn).reward_and_airpuff = 0; % reward-responsive AND airpuff-responsive
    agg.(rn).all_three        = 0;   % all three
    agg.(rn).none             = 0;   % none of the three
end

% Per-session proportion matrices for session-level statistics.
% Rows appended only for sessions where that region has >=1 neuron.
% enc columns      : [odor, reward, airpuff] / total neurons
% odor_out columns : [odor_reward, odor_airpuff, odor_neutral] / odor neurons
sess_props = struct();
for r = 1:length(regions)
    sess_props.(regions{r}).enc      = [];
    sess_props.(regions{r}).odor_out = [];
end



fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('Odor / Reward / Airpuff encoding — %% of all CA1 neurons by region\n');
fprintf('%s\n\n', repmat('=', 1, 78));

%% -------------------------------------------------------------------------
%% Main loop
%% -------------------------------------------------------------------------

for iFolder = 1:length(folder_list)

    cd(folder_list{iFolder});

    % Per-session CA1 threshold (mirrors v4)
    if strcmp(folder_list{iFolder}, 'D:\vCA1\M656\preprocessed\M656_2026_03_25')
        CA1_DEPTH_THRESH_UM = 400;
        fprintf('  [info] Increased CA1 depth threshold for M656_2026_03_25\n');
    else
        CA1_DEPTH_THRESH_UM = 300;
    end

    LoadExpKeys;

    % ---- Odor -> outcome map (mirrors swr_coactivity_same_vs_diff_v2.m) -----
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
    % STEP 1: Load clean_units and apply depth correction
    % (identical to v4)
    % ------------------------------------------------------------------
    record_depth    = ExpKeys.probe1_Depth;
    cfg_md.maxdepth = 1000 * record_depth;

    S = load_NPX_spike('clean_units_imec0.mat');
    S.usr.recording_depths = (cfg_md.maxdepth - S.usr.recording_depths) ...
                              .* cos(deg2rad(ExpKeys.probe1_insertion_roll));
    S.usr.probeID = repmat({'imec0'}, 1, length(S.label));

    if isfile('clean_units_imec1.mat')
        record_depth2 = ExpKeys.probe2_Depth;
        S2 = load_NPX_spike('clean_units_imec1.mat');
        S2.usr.recording_depths = (1000 * record_depth2 - S2.usr.recording_depths) ...
                                   .* cos(deg2rad(ExpKeys.probe2_insertion_roll));
        S2.usr.probeID = repmat({'imec1'}, 1, length(S2.label));

        S.label                = cat(2, S.label,                S2.label);
        S.t                    = cat(2, S.t,                    S2.t);
        S.usr.recording_depths = cat(2, S.usr.recording_depths, S2.usr.recording_depths);
        S.usr.shank_ids        = cat(2, S.usr.shank_ids,        S2.usr.shank_ids);
        S.usr.probeID          = cat(2, S.usr.probeID,          S2.usr.probeID);
        if isfield(S.usr, 'channel_ids') && isfield(S2.usr, 'channel_ids')
            S.usr.channel_ids  = cat(2, S.usr.channel_ids,      S2.usr.channel_ids);
        end
        clear S2;
    end

    nCells = length(S.label);

    % ------------------------------------------------------------------
    % STEP 2: Dorsal / ventral classification (mirrors v4 thresholds)
    % ------------------------------------------------------------------
    DH_threshold_imec0 = 3000;
    DH_threshold_imec1 = 3000;

    S.usr.HC_location = cell(1, nCells);
    for iCell = 1:nCells
        if strcmp(S.usr.probeID{iCell}, 'imec0')
            thresh = DH_threshold_imec0;
        else
            thresh = DH_threshold_imec1;
        end
        if S.usr.recording_depths(iCell) < thresh
            S.usr.HC_location{iCell} = 'Dorsal';
        else
            S.usr.HC_location{iCell} = 'Ventral';
        end
    end

    % ------------------------------------------------------------------
    % STEP 3: CA1 depth-band filter + decorrelation mask
    % Loop over SWR files exactly as in v4. For each file, find the best
    % SWR channel depth and mark neurons within ±CA1_DEPTH_THRESH_UM as
    % CA1-confirmed. A neuron is kept if it passes on any SWR file for
    % its probe. Decorr mask applied after CA1 filter, also as in v4.
    % ------------------------------------------------------------------
    ca1_confirmed = false(1, nCells);   % grows as SWR files are processed

    AllSWR = FindFiles('*SWR.mat');

    for currFile = 1:length(AllSWR)

        swrFile = AllSWR{currFile};

        % Identify probe (probeNum matches v4: 1=imec0, 2=imec1)
        if contains(swrFile, 'imec0')
            probeID  = 'imec0';
            probeNum = 1;
        elseif contains(swrFile, 'imec1')
            probeID  = 'imec1';
            probeNum = 2;
            if ~isfile('imec1_clean_lfp.mat')
                fprintf('  [skip] imec1 LFP not found for %s\n', swrFile);
                continue;
            end
        else
            fprintf('  [skip] Cannot determine probe from filename: %s\n', swrFile);
            continue;
        end

        % Identify region (lowercase to match v4 / selectSWRLfp)
        if contains(swrFile, 'dorsal')
            regionStr = 'dorsal';
        elseif contains(swrFile, 'ventral')
            regionStr = 'ventral';
        else
            fprintf('  [skip] Cannot determine dorsal/ventral from filename: %s\n', swrFile);
            continue;
        end

        % Check SWR count before loading LFP (skip if too few)
        tmp = load(swrFile);
        if strcmpi(regionStr, 'dorsal')
            tmpSWR = tmp.dorsal_swr;
        else
            tmpSWR = tmp.ventral_swr;
        end
        if length(tmpSWR.iv.tstart) < SWR_MIN_COUNT
            fprintf('  [skip] Too few SWRs (%d) in %s\n', ...
                    length(tmpSWR.iv.tstart), swrFile);
            continue;
        end
        clear tmp tmpSWR;

        % Load LFP and find best SWR channel depth (identical to v4)
        cfg_lfp        = [];
        cfg_lfp.probe  = probeNum;
        cfg_lfp.region = regionStr;

        csc        = loadLFP(cfg_lfp, cfg_md, ExpKeys);
        csc        = restrict(csc, 0, ExpKeys.sessend);
        lfp        = selectSWRLfp(csc, ExpKeys, cfg_lfp);
        swrChDepth = lfp.SWR.cfg.depths;

        fprintf('  [info] %s | best SWR ch depth = %.1f µm  (CA1 band: %.1f – %.1f µm)\n', ...
                swrFile, swrChDepth, ...
                swrChDepth - CA1_DEPTH_THRESH_UM, ...
                swrChDepth + CA1_DEPTH_THRESH_UM);

        % Mark neurons on this probe within the CA1 band
        probeIdx    = find(strcmp(S.usr.probeID, probeID));
        depthsProbe = S.usr.recording_depths(probeIdx);
        bandMask    = abs(depthsProbe - swrChDepth) <= CA1_DEPTH_THRESH_UM;
        ca1_confirmed(probeIdx(bandMask)) = true;

        fprintf('  [info] %s %s: %d / %d probe cells pass CA1 depth filter\n', ...
                probeID, regionStr, sum(bandMask), length(probeIdx));

    end  % SWR file loop

    % Apply CA1 filter
    ca1_idx                = find(ca1_confirmed);
    S.label                = S.label(ca1_idx);
    S.t                    = S.t(ca1_idx);
    S.usr.recording_depths = S.usr.recording_depths(ca1_idx);
    S.usr.shank_ids        = S.usr.shank_ids(ca1_idx);
    S.usr.probeID          = S.usr.probeID(ca1_idx);
    S.usr.HC_location      = S.usr.HC_location(ca1_idx);
    if isfield(S.usr, 'channel_ids')
        S.usr.channel_ids  = S.usr.channel_ids(ca1_idx);
    end
    nCells = length(S.label);

    fprintf('  [info] %d / %d total neurons pass CA1 depth filter\n', ...
            nCells, length(ca1_confirmed));

    if nCells == 0
        warning('  [warn] No CA1 neurons remain after CA1 filter in %s — skipping.\n', ...
                folder_list{iFolder});
        clearvars -except folder_list regions outcomes agg sess_props SWR_MIN_COUNT
        continue;
    end

    % ------------------------------------------------------------------
    % STEP 4: Decorrelation mask (applied after CA1 filter, as in v4)
    % ------------------------------------------------------------------
    decorr_keep = true(1, nCells);
    for iCell = 1:nCells
        probeID   = S.usr.probeID{iCell};
        shankID   = S.usr.shank_ids(iCell);
        fieldName = sprintf('%s_shank%d_decorr_times', probeID, shankID);
        if ~isfield(ExpKeys, fieldName), continue; end
        val = ExpKeys.(fieldName);
        if iscell(val), val = val{1}; end
        if isnumeric(val) && isscalar(val) && val >= 0
            decorr_keep(iCell) = false;
        end
    end

    n_removed = sum(~decorr_keep);
    if n_removed > 0
        fprintf('  [info] Decorr mask removed %d / %d CA1 neurons\n', n_removed, nCells);
    end

    keep_idx               = find(decorr_keep);
    S.label                = S.label(keep_idx);
    S.t                    = S.t(keep_idx);
    S.usr.recording_depths = S.usr.recording_depths(keep_idx);
    S.usr.shank_ids        = S.usr.shank_ids(keep_idx);
    S.usr.probeID          = S.usr.probeID(keep_idx);
    S.usr.HC_location      = S.usr.HC_location(keep_idx);
    nCells                 = length(S.label);

    if nCells == 0
        warning('  [warn] No neurons remain after decorr mask in %s — skipping.\n', ...
                folder_list{iFolder});
        clearvars -except folder_list regions outcomes agg sess_props SWR_MIN_COUNT
        continue;
    end

    % ------------------------------------------------------------------
    % STEP 5: Read neuron_summary.csv; collect encoding-positive ids
    % ------------------------------------------------------------------
    csv_path = 'neuron_summary.csv';
    if ~isfile(csv_path)
        warning('CSV not found, skipping: %s', folder_list{iFolder});
        clearvars -except folder_list regions outcomes agg sess_props SWR_MIN_COUNT
        continue;
    end
    data = readtable(csv_path, 'VariableNamingRule', 'preserve');

    % Build label map for the surviving CA1 neurons
    allUnitsMap = containers.Map(S.label, num2cell(1:nCells));

    odor_ids          = {};
    reward_ids        = {};
    airpuff_ids       = {};
    odor_reward_ids   = {};   % odor-selective AND reward outcome
    odor_airpuff_ids  = {};   % odor-selective AND airpuff outcome
    odor_neutral_ids  = {};   % odor-selective AND neutral outcome

    for row = 1:height(data)

        global_id_raw   = data.global_id{row};
        global_id_split = split(global_id_raw, '.');
        global_id       = strcat(global_id_split{1}, '_', global_id_split{3});

        % Skip neurons not in the surviving CA1 set
        if ~isKey(allUnitsMap, global_id)
            continue;
        end

        region  = lower(strtrim(data.region{row}));
        skip_fr = data.skip_low_baseline_fr{row};

        if ismember(region, {'outside', 'outside_ca1'})
            continue;
        end
        if strcmpi(skip_fr, 'Yes')
            continue;
        end

        odor_sel      = data.odor_selectivity{row};
        odor_resp_raw = data.odor_responsivity{row};
        reward_resp   = data.reward_responsive{row};
        airpuff_resp  = data.airpuff_responsive{row};

        odor_sel_is_empty = isempty(odor_sel) || strcmpi(odor_sel, 'None') || ...
            (isscalar(odor_sel) && isnumeric(odor_sel) && isnan(odor_sel));
        odor_resp_is_empty = isempty(odor_resp_raw) || strcmpi(odor_resp_raw, 'None') || ...
            (isscalar(odor_resp_raw) && isnumeric(odor_resp_raw) && isnan(odor_resp_raw));

        is_odor_selective = ~odor_sel_is_empty;
        if ~is_odor_selective && ~odor_resp_is_empty
            responsive_odors = strtrim(strsplit(odor_resp_raw, ','));
            responsive_odors = responsive_odors(~cellfun('isempty', responsive_odors));
            if isscalar(responsive_odors)
                is_odor_selective = true;
            end
        end

        is_reward  = strcmpi(reward_resp,  'Yes');
        is_airpuff = strcmpi(airpuff_resp, 'Yes');

        if is_odor_selective,  odor_ids{end+1}    = global_id; end %#ok<AGROW>
        if is_reward,          reward_ids{end+1}   = global_id; end %#ok<AGROW>
        if is_airpuff,         airpuff_ids{end+1}  = global_id; end %#ok<AGROW>

        % Odor-selective neurons broken down by outcome category.
        % Look up the selective odor(s) in the ExpKeys-derived map.
        % Priority when a neuron is selective for odors of mixed outcomes
        % (shouldn't happen with single-odor promotion, but kept for safety):
        %   Reward > Airpuff > Neutral
        if is_odor_selective

            % Collect the odor name(s) this neuron is selective for.
            % odor_selectivity holds the promoted single odor; if empty,
            % fall back to odor_responsivity (single-responsive promotion).
            if ~odor_sel_is_empty
                candidate_odors = strtrim(strsplit(odor_sel, ','));
            else
                candidate_odors = strtrim(strsplit(odor_resp_raw, ','));
            end
            candidate_odors = candidate_odors(~cellfun('isempty', candidate_odors));

            % Map each odor to its outcome; resolve with priority order
            og = 'Neutral';   % default if no odor found in map
            for iOdor = 1:length(candidate_odors)
                o = candidate_odors{iOdor};
                if isKey(odor_to_outcome, o)
                    mapped = odor_to_outcome(o);
                    % Upgrade priority: Reward > Airpuff > Neutral
                    if strcmp(mapped, 'Reward')
                        og = 'Reward';
                        break;   % can't do better
                    elseif strcmp(mapped, 'Airpuff') && ~strcmp(og, 'Reward')
                        og = 'Airpuff';
                    end
                    % Neutral: leave og as-is (already Neutral or higher)
                end
            end

            switch lower(og)
                case 'reward'
                    odor_reward_ids{end+1}  = global_id; %#ok<AGROW>
                case 'airpuff'
                    odor_airpuff_ids{end+1} = global_id; %#ok<AGROW>
                otherwise
                    odor_neutral_ids{end+1} = global_id; %#ok<AGROW>
            end
        end

    end  % row loop

    % ------------------------------------------------------------------
    % STEP 6: Count totals and encoding-positive cells per region
    % ------------------------------------------------------------------
    sess = struct();
    for r = 1:length(regions)
        rn = regions{r};
        sess.(rn).total        = 0;
        sess.(rn).odor         = 0;
        sess.(rn).reward       = 0;
        sess.(rn).airpuff      = 0;
        sess.(rn).odor_reward  = 0;
        sess.(rn).odor_airpuff = 0;
        sess.(rn).odor_neutral = 0;
        % Overlap counts
        sess.(rn).odor_only         = 0;
        sess.(rn).reward_only       = 0;
        sess.(rn).airpuff_only      = 0;
        sess.(rn).odor_and_reward   = 0;
        sess.(rn).odor_and_airpuff  = 0;
        sess.(rn).reward_and_airpuff = 0;
        sess.(rn).all_three         = 0;
        sess.(rn).none              = 0;
    end

    for iCell = 1:nCells
        lbl = S.label{iCell};
        loc = S.usr.HC_location{iCell};
        if ~ismember(loc, regions), continue; end

        sess.(loc).total = sess.(loc).total + 1;

        if ismember(lbl, odor_ids),    sess.(loc).odor    = sess.(loc).odor    + 1; end
        if ismember(lbl, reward_ids),  sess.(loc).reward  = sess.(loc).reward  + 1; end
        if ismember(lbl, airpuff_ids), sess.(loc).airpuff = sess.(loc).airpuff + 1; end

        if ismember(lbl, odor_reward_ids),  sess.(loc).odor_reward  = sess.(loc).odor_reward  + 1; end
        if ismember(lbl, odor_airpuff_ids), sess.(loc).odor_airpuff = sess.(loc).odor_airpuff + 1; end
        if ismember(lbl, odor_neutral_ids), sess.(loc).odor_neutral = sess.(loc).odor_neutral + 1; end

        % Overlap bins — classify each neuron into one of 8 exclusive groups
        has_o = ismember(lbl, odor_ids);
        has_r = ismember(lbl, reward_ids);
        has_a = ismember(lbl, airpuff_ids);
        if     has_o && has_r && has_a
            sess.(loc).all_three          = sess.(loc).all_three          + 1;
        elseif has_o && has_r
            sess.(loc).odor_and_reward    = sess.(loc).odor_and_reward    + 1;
        elseif has_o && has_a
            sess.(loc).odor_and_airpuff   = sess.(loc).odor_and_airpuff   + 1;
        elseif has_r && has_a
            sess.(loc).reward_and_airpuff = sess.(loc).reward_and_airpuff + 1;
        elseif has_o
            sess.(loc).odor_only          = sess.(loc).odor_only          + 1;
        elseif has_r
            sess.(loc).reward_only        = sess.(loc).reward_only        + 1;
        elseif has_a
            sess.(loc).airpuff_only       = sess.(loc).airpuff_only       + 1;
        else
            sess.(loc).none               = sess.(loc).none               + 1;
        end
    end

    % ------------------------------------------------------------------
    % Print session summary — encoding overview
    % ------------------------------------------------------------------
    [~, sess_name] = fileparts(folder_list{iFolder});
    fprintf('Session: %s\n', sess_name);
    fprintf('%-8s  %5s  %16s  %18s  %18s\n', ...
        'Region', 'N', 'Odor-selective', 'Reward-responsive', 'Airpuff-responsive');
    fprintf('%s\n', repmat('-', 1, 72));

    for r = 1:length(regions)
        rn = regions{r};
        N  = sess.(rn).total;
        if N == 0
            fprintf('%-8s  %5d  (no neurons)\n', rn, 0);
        else
            fprintf('%-8s  %5d  %4d (%5.1f%%)   %4d (%5.1f%%)     %4d (%5.1f%%)\n', ...
                rn, N, ...
                sess.(rn).odor,    100 * sess.(rn).odor    / N, ...
                sess.(rn).reward,  100 * sess.(rn).reward  / N, ...
                sess.(rn).airpuff, 100 * sess.(rn).airpuff / N);
        end
    end

    % ------------------------------------------------------------------
    % Print session outcome-stratified odor breakdown
    % ------------------------------------------------------------------
    fprintf('\n  Odor-selective neurons by outcome category:\n');
    fprintf('  %-8s  %5s  %5s (%%odor)  %5s (%%odor)  %5s (%%odor)\n', ...
        'Region', 'Odor', 'Reward', 'Airpuff', 'Neutral');
    fprintf('  %s\n', repmat('-', 1, 62));
    for r = 1:length(regions)
        rn  = regions{r};
        nOd = sess.(rn).odor;
        if nOd == 0
            fprintf('  %-8s  %5d  (no odor-selective neurons)\n', rn, 0);
        else
            fprintf('  %-8s  %5d  %5d (%5.1f%%)  %5d (%5.1f%%)  %5d (%5.1f%%)\n', ...
                rn, nOd, ...
                sess.(rn).odor_reward,  100 * sess.(rn).odor_reward  / nOd, ...
                sess.(rn).odor_airpuff, 100 * sess.(rn).odor_airpuff / nOd, ...
                sess.(rn).odor_neutral, 100 * sess.(rn).odor_neutral / nOd);
        end
    end
    fprintf('\n');

    % ------------------------------------------------------------------
    % Print session overlap summary
    % ------------------------------------------------------------------
    fprintf('\n  Encoding overlap (exclusive bins, %% of total CA1 neurons):\n');
    fprintf('  %-8s  %5s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s\n', ...
        'Region', 'N', 'Odor only', 'Reward only', 'Airpuff only', ...
        'Odor+Reward', 'Odor+Airpuff', 'Rew+Airpuff', 'All three', 'None');
    fprintf('  %s\n', repmat('-', 1, 122));
    for r = 1:length(regions)
        rn = regions{r};
        N  = sess.(rn).total;
        if N == 0
            fprintf('  %-8s  %5d  (no neurons)\n', rn, 0);
        else
            pct = @(x) 100 * x / N;
            fprintf('  %-8s  %5d  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)\n', ...
                rn, N, ...
                sess.(rn).odor_only,          pct(sess.(rn).odor_only), ...
                sess.(rn).reward_only,        pct(sess.(rn).reward_only), ...
                sess.(rn).airpuff_only,       pct(sess.(rn).airpuff_only), ...
                sess.(rn).odor_and_reward,    pct(sess.(rn).odor_and_reward), ...
                sess.(rn).odor_and_airpuff,   pct(sess.(rn).odor_and_airpuff), ...
                sess.(rn).reward_and_airpuff, pct(sess.(rn).reward_and_airpuff), ...
                sess.(rn).all_three,          pct(sess.(rn).all_three), ...
                sess.(rn).none,               pct(sess.(rn).none));
        end
    end
    fprintf('\n');

    % ------------------------------------------------------------------
    % Accumulate aggregate totals and per-session proportions
    % ------------------------------------------------------------------
    for r = 1:length(regions)
        rn = regions{r};
        agg.(rn).total        = agg.(rn).total        + sess.(rn).total;
        agg.(rn).odor         = agg.(rn).odor         + sess.(rn).odor;
        agg.(rn).reward       = agg.(rn).reward       + sess.(rn).reward;
        agg.(rn).airpuff      = agg.(rn).airpuff      + sess.(rn).airpuff;
        agg.(rn).odor_reward  = agg.(rn).odor_reward  + sess.(rn).odor_reward;
        agg.(rn).odor_airpuff = agg.(rn).odor_airpuff + sess.(rn).odor_airpuff;
        agg.(rn).odor_neutral = agg.(rn).odor_neutral + sess.(rn).odor_neutral;
        % Overlap
        agg.(rn).odor_only          = agg.(rn).odor_only          + sess.(rn).odor_only;
        agg.(rn).reward_only        = agg.(rn).reward_only        + sess.(rn).reward_only;
        agg.(rn).airpuff_only       = agg.(rn).airpuff_only       + sess.(rn).airpuff_only;
        agg.(rn).odor_and_reward    = agg.(rn).odor_and_reward    + sess.(rn).odor_and_reward;
        agg.(rn).odor_and_airpuff   = agg.(rn).odor_and_airpuff   + sess.(rn).odor_and_airpuff;
        agg.(rn).reward_and_airpuff = agg.(rn).reward_and_airpuff + sess.(rn).reward_and_airpuff;
        agg.(rn).all_three          = agg.(rn).all_three          + sess.(rn).all_three;
        agg.(rn).none               = agg.(rn).none               + sess.(rn).none;
    end

    % Append session-level proportions for each region that has neurons
    for r = 1:length(regions)
        rn = regions{r};
        N   = sess.(rn).total;
        nOd = sess.(rn).odor;
        if N > 0
            sess_props.(rn).enc(end+1, :) = ...
                [sess.(rn).odor / N, sess.(rn).reward / N, sess.(rn).airpuff / N];
        end
        if nOd > 0
            sess_props.(rn).odor_out(end+1, :) = ...
                [sess.(rn).odor_reward  / nOd, ...
                 sess.(rn).odor_airpuff / nOd, ...
                 sess.(rn).odor_neutral / nOd];
        end
    end

    clearvars -except folder_list regions outcomes agg sess_props SWR_MIN_COUNT
end  % folder loop

%% -------------------------------------------------------------------------
%% Aggregate summary — encoding overview
%% -------------------------------------------------------------------------

fprintf('%s\n', repmat('=', 1, 72));
fprintf('AGGREGATE — all sessions\n');
fprintf('%-8s  %5s  %16s  %18s  %18s\n', ...
    'Region', 'N', 'Odor-selective', 'Reward-responsive', 'Airpuff-responsive');
fprintf('%s\n', repmat('-', 1, 72));

for r = 1:length(regions)
    rn = regions{r};
    N  = agg.(rn).total;
    if N == 0
        fprintf('%-8s  %5d  (no neurons)\n', rn, 0);
    else
        fprintf('%-8s  %5d  %4d (%5.1f%%)   %4d (%5.1f%%)     %4d (%5.1f%%)\n', ...
            rn, N, ...
            agg.(rn).odor,    100 * agg.(rn).odor    / N, ...
            agg.(rn).reward,  100 * agg.(rn).reward  / N, ...
            agg.(rn).airpuff, 100 * agg.(rn).airpuff / N);
    end
end
fprintf('%s\n\n', repmat('=', 1, 72));

%% -------------------------------------------------------------------------
%% Aggregate — odor-selective neurons by outcome category
%% -------------------------------------------------------------------------

fprintf('%s\n', repmat('=', 1, 72));
fprintf('AGGREGATE — odor-selective neurons by outcome category\n');
fprintf('%-8s  %5s  %5s (%%odor)  %5s (%%odor)  %5s (%%odor)\n', ...
    'Region', 'Odor', 'Reward', 'Airpuff', 'Neutral');
fprintf('%s\n', repmat('-', 1, 62));

for r = 1:length(regions)
    rn  = regions{r};
    nOd = agg.(rn).odor;
    if nOd == 0
        fprintf('%-8s  %5d  (no odor-selective neurons)\n', rn, 0);
    else
        fprintf('%-8s  %5d  %5d (%5.1f%%)  %5d (%5.1f%%)  %5d (%5.1f%%)\n', ...
            rn, nOd, ...
            agg.(rn).odor_reward,  100 * agg.(rn).odor_reward  / nOd, ...
            agg.(rn).odor_airpuff, 100 * agg.(rn).odor_airpuff / nOd, ...
            agg.(rn).odor_neutral, 100 * agg.(rn).odor_neutral / nOd);
    end
end
fprintf('%s\n\n', repmat('=', 1, 72));

%% -------------------------------------------------------------------------
%% Aggregate — encoding overlap (exclusive bins)
%% -------------------------------------------------------------------------

fprintf('%s\n', repmat('=', 1, 72));
fprintf('AGGREGATE — encoding overlap (exclusive bins, %% of total CA1 neurons)\n');
fprintf('%-8s  %5s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s\n', ...
    'Region', 'N', 'Odor only', 'Reward only', 'Airpuff only', ...
    'Odor+Reward', 'Odor+Airpuff', 'Rew+Airpuff', 'All three', 'None');
fprintf('%s\n', repmat('-', 1, 122));

for r = 1:length(regions)
    rn = regions{r};
    N  = agg.(rn).total;
    if N == 0
        fprintf('%-8s  %5d  (no neurons)\n', rn, 0);
    else
        pct = @(x) 100 * x / N;
        fprintf('%-8s  %5d  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)  %4d (%5.1f%%)\n', ...
            rn, N, ...
            agg.(rn).odor_only,          pct(agg.(rn).odor_only), ...
            agg.(rn).reward_only,        pct(agg.(rn).reward_only), ...
            agg.(rn).airpuff_only,       pct(agg.(rn).airpuff_only), ...
            agg.(rn).odor_and_reward,    pct(agg.(rn).odor_and_reward), ...
            agg.(rn).odor_and_airpuff,   pct(agg.(rn).odor_and_airpuff), ...
            agg.(rn).reward_and_airpuff, pct(agg.(rn).reward_and_airpuff), ...
            agg.(rn).all_three,          pct(agg.(rn).all_three), ...
            agg.(rn).none,               pct(agg.(rn).none));
    end
end
fprintf('%s\n\n', repmat('=', 1, 72));

%% -------------------------------------------------------------------------
%% Statistical tests — Dorsal vs. Ventral
%%
%% Three complementary analyses:
%%
%%  A) Neuron-level (aggregate): Fisher's exact on pooled 2x2 counts.
%%     Treats every neuron as an independent observation.
%%     Bonferroni for 3 encoding types: alpha = 0.017.
%%
%%  B) Session-level: Wilcoxon rank-sum (unpaired) on per-session proportions.
%%     Sessions contribute one proportion each; dorsal and ventral sessions
%%     are almost never the same recording so this is unpaired.
%%     Bonferroni for 3 encoding types: alpha = 0.017.
%%
%%  C) Odor-outcome category: do dorsal and ventral differ in the
%%     distribution of odor-selective neurons across Reward/Airpuff/Neutral?
%%     - Chi-square test on the 2x3 contingency table (omnibus).
%%     - Pairwise Fisher's exact for each outcome column vs the other two
%%       (3 pairs, Bonferroni alpha = 0.017).
%%     Session-level Wilcoxon rank-sum on per-session odor_out proportions
%%     (same 3 pairwise comparisons).
%% -------------------------------------------------------------------------

ALPHA      = 0.05;
ALPHA_BONF = ALPHA / 3;   % 3 encoding types or 3 outcome pairs

fprintf('%s\n', repmat('=', 1, 72));
fprintf('STATISTICS — Dorsal vs. Ventral\n');
fprintf('%s\n\n', repmat('=', 1, 72));

% ---- helpers ----
enc_labels    = {'Odor-selective',  'Reward-responsive', 'Airpuff-responsive'};
enc_fields    = {'odor',            'reward',             'airpuff'};
out_labels    = {'Reward-odor',     'Airpuff-odor',       'Neutral-odor'};
out_agg_flds  = {'odor_reward',     'odor_airpuff',       'odor_neutral'};

props_D = sess_props.Dorsal;
props_V = sess_props.Ventral;

% ==========================================================================
% A) Neuron-level Fisher's exact — encoding proportions
% ==========================================================================
fprintf('A) NEURON-LEVEL — Fisher''s exact (aggregate counts)\n');
fprintf('   Bonferroni-corrected alpha = %.3f for %d tests\n\n', ALPHA_BONF, 3);
fprintf('   %-20s  %8s  %8s  %10s\n', '', 'Dorsal', 'Ventral', 'p');
fprintf('   %s\n', repmat('-', 1, 52));

for k = 1:3
    nd_enc = agg.Dorsal.(enc_fields{k});
    nd_tot = agg.Dorsal.total;
    nv_enc = agg.Ventral.(enc_fields{k});
    nv_tot = agg.Ventral.total;

    pct_d = 100 * nd_enc / nd_tot;
    pct_v = 100 * nv_enc / nv_tot;

    [~, p] = fishertest([nd_enc, nd_tot - nd_enc; nv_enc, nv_tot - nv_enc]);
    if p < ALPHA_BONF, sig = '*'; else, sig = 'ns'; end

    fprintf('   %-20s  %5.1f%%       %5.1f%%    p = %.4f  [%s]\n', ...
            enc_labels{k}, pct_d, pct_v, p, sig);
end
fprintf('\n');

% ==========================================================================
% B) Session-level Wilcoxon rank-sum — encoding proportions
% ==========================================================================
fprintf('B) SESSION-LEVEL — Wilcoxon rank-sum on per-session proportions\n');
fprintf('   (unpaired; sessions rarely shared across regions)\n');
fprintf('   Bonferroni-corrected alpha = %.3f for %d tests\n\n', ALPHA_BONF, 3);
fprintf('   %-20s  %16s  %16s  %10s\n', '', ...
        'Dorsal (n sess)', 'Ventral (n sess)', 'p');
fprintf('   %s\n', repmat('-', 1, 72));

for k = 1:3
    d_vec = props_D.enc(:, k);
    v_vec = props_V.enc(:, k);

    fprintf('   %-20s  %5.1f%% (n=%d)     %5.1f%% (n=%d)  ', ...
            enc_labels{k}, ...
            mean(d_vec)*100, numel(d_vec), ...
            mean(v_vec)*100, numel(v_vec));

    if numel(d_vec) < 2 || numel(v_vec) < 2
        fprintf('  [skip — need >=2 sessions per region]\n');
    else
        [p, ~] = ranksum(d_vec, v_vec);
        if p < ALPHA_BONF, sig = '*'; else, sig = 'ns'; end
        fprintf('p = %.4f  [%s]\n', p, sig);
    end
end
fprintf('\n');

% ==========================================================================
% C) Odor-outcome category — dorsal vs. ventral
% ==========================================================================
fprintf('C) ODOR-OUTCOME CATEGORY — Dorsal vs. Ventral\n');
fprintf('   Among odor-selective neurons: Reward / Airpuff / Neutral outcome\n\n');

% --- C1: omnibus chi-square on 2x3 table (neuron-level) ---
fprintf('   C1) Omnibus chi-square on aggregate 2x3 contingency table\n');
ct2x3 = [ agg.Dorsal.odor_reward,  agg.Dorsal.odor_airpuff,  agg.Dorsal.odor_neutral; ...
          agg.Ventral.odor_reward, agg.Ventral.odor_airpuff, agg.Ventral.odor_neutral];
fprintf('        Dorsal  (n=%d odor): %d reward, %d airpuff, %d neutral\n', ...
        agg.Dorsal.odor,  ct2x3(1,1), ct2x3(1,2), ct2x3(1,3));
fprintf('        Ventral (n=%d odor): %d reward, %d airpuff, %d neutral\n', ...
        agg.Ventral.odor, ct2x3(2,1), ct2x3(2,2), ct2x3(2,3));

% Inline chi-square: expected counts, stat, df, p
row_sums  = sum(ct2x3, 2);
col_sums  = sum(ct2x3, 1);
n_total   = sum(col_sums);
expected  = (row_sums * col_sums) / n_total;
chi2stat  = sum(sum((ct2x3 - expected).^2 ./ expected));
df_chi    = (size(ct2x3,1)-1) * (size(ct2x3,2)-1);
p_chi     = 1 - chi2cdf(chi2stat, df_chi);
fprintf('        Chi-square(%d) = %.3f,  p = %.4f\n\n', df_chi, chi2stat, p_chi);

% --- C2: pairwise Fisher's exact per outcome column (neuron-level) ---
fprintf('   C2) Pairwise Fisher''s exact per outcome (neuron-level)\n');
fprintf('       Bonferroni-corrected alpha = %.3f for %d tests\n\n', ALPHA_BONF, 3);
fprintf('       %-16s  %8s  %8s  %10s\n', 'Outcome', 'Dorsal', 'Ventral', 'p');
fprintf('       %s\n', repmat('-', 1, 48));

for k = 1:3
    nd_enc = agg.Dorsal.(out_agg_flds{k});
    nd_tot = agg.Dorsal.odor;
    nv_enc = agg.Ventral.(out_agg_flds{k});
    nv_tot = agg.Ventral.odor;

    if nd_tot == 0 || nv_tot == 0
        fprintf('       %-16s  (insufficient data)\n', out_labels{k});
        continue;
    end

    pct_d = 100 * nd_enc / nd_tot;
    pct_v = 100 * nv_enc / nv_tot;
    [~, p] = fishertest([nd_enc, nd_tot - nd_enc; nv_enc, nv_tot - nv_enc]);
    if p < ALPHA_BONF, sig = '*'; else, sig = 'ns'; end

    fprintf('       %-16s  %5.1f%%       %5.1f%%    p = %.4f  [%s]\n', ...
            out_labels{k}, pct_d, pct_v, p, sig);
end
fprintf('\n');

% --- C3: session-level Wilcoxon rank-sum per outcome proportion ---
fprintf('   C3) Session-level Wilcoxon rank-sum per outcome proportion\n');
fprintf('       (unpaired; proportion of odor-selective neurons per session)\n');
fprintf('       Bonferroni-corrected alpha = %.3f for %d tests\n\n', ALPHA_BONF, 3);
fprintf('       %-16s  %16s  %16s  %10s\n', 'Outcome', ...
        'Dorsal (n sess)', 'Ventral (n sess)', 'p');
fprintf('       %s\n', repmat('-', 1, 68));

for k = 1:3
    d_vec = props_D.odor_out(:, k);
    v_vec = props_V.odor_out(:, k);

    if isempty(d_vec) || isempty(v_vec)
        fprintf('       %-16s  (no sessions with odor-selective neurons in one region)\n', ...
                out_labels{k});
        continue;
    end

    fprintf('       %-16s  %5.1f%% (n=%d)     %5.1f%% (n=%d)  ', ...
            out_labels{k}, ...
            mean(d_vec)*100, numel(d_vec), ...
            mean(v_vec)*100, numel(v_vec));

    if numel(d_vec) < 2 || numel(v_vec) < 2
        fprintf('  [skip — need >=2 sessions per region]\n');
    else
        [p, ~] = ranksum(d_vec, v_vec);
        if p < ALPHA_BONF, sig = '*'; else, sig = 'ns'; end
        fprintf('p = %.4f  [%s]\n', p, sig);
    end
end
fprintf('\n');

fprintf('%s\n\n', repmat('=', 1, 72));
fprintf('Script completed.\n');

end  % main function

%% =========================================================================
%% Local functions (copied from swr_participation_stimulus_encoding_v4.m)
%% =========================================================================

function lfp = selectSWRLfp(csc, ExpKeys, cfg)
% selectSWRLfp  Extract the SWR channel for the chosen probe region.
    if strcmpi(cfg.region, 'dorsal')
        swrCh = ExpKeys.imec0_dorsal_SWR_channel;
    else
        swrCh = ExpKeys.imec0_ventral_SWR_channel;
    end
    swrIdx             = find(cellfun(@(x) strcmp(x, swrCh), csc.label), 1);
    lfp.SWR            = csc;
    lfp.SWR.data       = lfp.SWR.data(swrIdx, :);
    lfp.SWR.cfg.depths = lfp.SWR.cfg.depths(swrIdx);
    lfp.SWR.label      = lfp.SWR.label(swrIdx);
end

function csc = loadLFP(cfg, cfg_md, ExpKeys)
% loadLFP  Load LFP for the requested probe, apply true-depth correction,
%          and sort channels dorsal (shallow) -> ventral (deep).
%
%   cfg.probe  : 1 = imec0, 2 = imec1
%   cfg_md     : struct with field maxdepth (µm)

    csc = tsd;

    if cfg.probe == 1
        load('imec0_clean_lfp.mat');
        csc.tvec          = imec0.lfp_tvec';
        csc.data          = double(imec0.lfp_traces');
        csc.label         = cellstr(string(imec0.channel_ids).strip().extractAfter('ap#')');
        csc.cfg.depths    = imec0.depths;
        csc.cfg.shank_ids = imec0.shank_ids;
        Fs = struct; Fs.SamplingFrequency = double(imec0.lfp_fs);
        csc.cfg.hdr{1,1}  = Fs;
        csc.cfg.depths    = (cfg_md.maxdepth - csc.cfg.depths) ...
                             .* cos(deg2rad(ExpKeys.probe1_insertion_roll));
    else
        load('imec1_clean_lfp.mat');
        csc.tvec          = imec1.lfp_tvec';
        csc.data          = double(imec1.lfp_traces');
        csc.label         = cellstr(string(imec1.channel_ids).strip().extractAfter('ap#')');
        csc.cfg.depths    = imec1.depths;
        csc.cfg.shank_ids = imec1.shank_ids;
        Fs = struct; Fs.SamplingFrequency = double(imec1.lfp_fs);
        csc.cfg.hdr{1,1}  = Fs;
        maxdepth2         = 1000 * ExpKeys.probe2_Depth;
        csc.cfg.depths    = (maxdepth2 - csc.cfg.depths) ...
                             .* cos(deg2rad(ExpKeys.probe2_insertion_roll));
    end

    [~, sort_idx]     = sort(csc.cfg.depths);
    csc.data          = csc.data(sort_idx, :);
    csc.label         = csc.label(sort_idx);
    csc.cfg.depths    = csc.cfg.depths(sort_idx);
    csc.cfg.shank_ids = csc.cfg.shank_ids(sort_idx);
end