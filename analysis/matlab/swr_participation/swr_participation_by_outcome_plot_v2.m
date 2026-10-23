%% Sharp Wave Ripple Participation - By Outcome / Trial Type Group
% Extension of swr_participation_stimulus_encoding_v4.m
%
% Tests whether neurons encoding different outcome categories (Reward /
% Airpuff / Neutral / Unclassified) participate in SWRs at different rates,
% and whether reward-encoding neurons differ from airpuff-encoding neurons.
%
% -------------------------------------------------------------------------
% NEURON CLASSIFICATION  (mirrors swr_coactivity_same_vs_diff_v2.m exactly)
%
%   Inclusion criterion (same as coactivity script):
%     At least one EXCITATORY / INCREASE response:
%       odor_dir    == 'excitatory' or 'increase'
%       reward_dir  == 'increase'
%       airpuff_dir == 'increase'
%
%   Trial-type set per neuron (union of):
%     - Outcome linked to its preferred odor via odor_to_outcome map
%       (Reward / Airpuff / Neutral), if odor_dir is defined
%     - 'Reward'  if reward_dir  ~= 'None'
%     - 'Airpuff' if airpuff_dir ~= 'None'
%
%   Primary outcome group (for pairwise comparisons):
%     Assigned by priority:  Reward > Neutral > Airpuff > Unclassified
%     A neuron whose trial-type set contains 'Reward' is labelled 'Reward'
%     even if it also encodes 'Airpuff'.
%
%   TrialTypeSet column stores the full set as a comma-separated string
%   (e.g., 'Reward,Airpuff') for downstream inspection.
%
% -------------------------------------------------------------------------
% ANALYSES
%
%   1. Outcome group (encoding neurons only):
%      Kruskal-Wallis across Reward / Airpuff / Neutral (3 groups).
%      Unclassified neurons are excluded to keep the omnibus test and
%      Bonferroni correction focused on the outcome-identity question.
%      Post-hoc: Wilcoxon rank-sum for every pair, Bonferroni-corrected
%      (at most 3 pairs -> alpha = 0.017).
%      This version subs a wilcoxn-rank sum for a Post-hoc Tukey-Kramer on ranks
%
%   2. Reward vs Airpuff (encoding neurons only):
%      Wilcoxon rank-sum, overall and per region.
%      Membership defined by trial-type SET containing the relevant
%      outcome; neurons encoding both are counted on both sides and flagged.
%
%   3. Encoding vs Unclassified (full population):
%      Wilcoxon rank-sum reprising swr_participation_stimulus_encoding_v4.m
%      so the broader encoding/non-encoding contrast is preserved.
%
% -------------------------------------------------------------------------
% DEPENDENCIES  (same as v4 participation script)
%   LoadExpKeys, FindFiles, load_NPX_spike, ResizeIV, restrict, iv,
%   SelectTS, selectSWRLfp (or equivalent), loadLFP (local function below)
%
% OUTPUT
%   celltable  - one row per qualifying neuron, extended with:
%                OutcomeGroup, TrialTypeSet, OdorOutcome,
%                RewardDir, AirpuffDir
%   Saved as *_SWR_participation_by_outcome.mat / .csv
%
% EG - June 2026
% v1: Extended from swr_participation_stimulus_encoding_v4.m
%     Trial-type logic from swr_coactivity_same_vs_diff_v2.m
% v2: Restructured statistics - Analysis 1 restricted to encoding neurons
%     (3-group KW); Analysis 3 added for encoding vs unclassified contrast

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

SWR_MIN_COUNT       = 20;    % min SWRs required to include a session/file
CA1_DEPTH_THRESH_UM = 300;   % ±µm around best SWR channel to define CA1 band

MIN_BASELINE_FR = 0;   % Hz; neurons below this are excluded

spikeFileBase = 'selective_permissive';

cfg_iv   = [];
cfg_iv.d = [0, 0];   % expand each SWR interval (s) on each side

%% =========================================================================
%% Initialize accumulators
%% =========================================================================

all_session           = {};
all_label             = {};
all_location          = {};
all_cell_type         = {};    % Odor_Selective / Outcome_Responsive / Odor_and_Outcome / Unclassified
all_outcome_group     = {};    % Primary outcome: Reward / Airpuff / Neutral / Unclassified
all_trial_type_set    = {};    % Full trial-type set as comma-separated string
all_odor_outcome      = {};    % Outcome linked to preferred odor (or 'None')
all_reward_dir        = {};
all_airpuff_dir       = {};
all_stimulus_encoding = [];    % logical: has ≥1 excitatory/increase response
all_pct_participation = [];
all_swr_count         = [];
all_ExpKeys           = {};

%% =========================================================================
%% Main loop - folders
%% =========================================================================

for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});

    % Session-specific CA1 depth threshold
    if strcmp(folder_list{iFolder,1}, 'D:\vCA1\M656\preprocessed\M656_2026_03_25')
        CA1_DEPTH_THRESH_UM = 400;
        fprintf('Increased CA1 depth threshold for M656_2026_03_25\n');
    else
        CA1_DEPTH_THRESH_UM = 300;
    end

    rng(2024);
    LoadExpKeys;
    all_ExpKeys{end+1} = ExpKeys;

    % ---- Depth calibration (probe 1) ----------------------------------------
    record_depth    = ExpKeys.probe1_Depth;
    cfg_md.maxdepth = 1000 * record_depth;

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

    % ---- Load all clean units -----------------------------------------------
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
        S.usr.mean_waveforms   = cat(2, S.usr.mean_waveforms,   S2.usr.mean_waveforms);
        S.usr.probeID          = cat(2, S.usr.probeID,          S2.usr.probeID);
        if isfield(S.usr, 'channel_ids') && isfield(S2.usr, 'channel_ids')
            S.usr.channel_ids  = cat(2, S.usr.channel_ids,      S2.usr.channel_ids);
        end
        clear S2;
    end

    % ---- Initialise selectivity fields to 'None' ----------------------------
    nCells = length(S.label);
    S.usr.odor        = repmat({'None'}, 1, nCells);
    S.usr.odor_dir    = repmat({'None'}, 1, nCells);
    S.usr.reward_dir  = repmat({'None'}, 1, nCells);
    S.usr.airpuff_dir = repmat({'None'}, 1, nCells);

    allUnitsMap = containers.Map(S.label, num2cell(1:nCells));

    selectiveFiles = { ...
        strcat(spikeFileBase, '_imec0_all.mat'), ...
        strcat(spikeFileBase, '_imec1_all.mat')  ...
        };
    for iSF = 1:length(selectiveFiles)
        if exist(selectiveFiles{iSF}, 'file') ~= 2
            fprintf('  [info] %s not found - probe cells treated as Unclassified\n', ...
                selectiveFiles{iSF});
            continue;
        end
        load(selectiveFiles{iSF});   % loads S_out
        for iSel = 1:length(S_out.label)
            lbl = S_out.label{iSel};
            if isKey(allUnitsMap, lbl)
                idx = allUnitsMap(lbl);
                S.usr.odor{idx}        = S_out.usr.odor{iSel};
                S.usr.odor_dir{idx}    = S_out.usr.odor_dir{iSel};
                S.usr.reward_dir{idx}  = S_out.usr.reward_dir{iSel};
                S.usr.airpuff_dir{idx} = S_out.usr.airpuff_dir{iSel};
            else
                fprintf('  [warn] Label %s in %s not found in clean_units - skipping\n', ...
                    lbl, selectiveFiles{iSF});
            end
        end
        clear S_out;
    end

    % ---- Dorsal / ventral classification ------------------------------------
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

    % ---- SWR file loop ------------------------------------------------------
    AllSWR    = FindFiles('*SWR.mat');
    sessionID = strcat(ExpKeys.subject, '-', ExpKeys.date);

    for currFile = 1:length(AllSWR)

        swrFile = AllSWR{currFile};

        % Identify probe
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

        % Identify region
        if contains(swrFile, 'dorsal')
            regionStr    = 'dorsal';
            swr_location = 'Dorsal';
        elseif contains(swrFile, 'ventral')
            regionStr    = 'ventral';
            swr_location = 'Ventral';
        else
            fprintf('  [skip] Cannot determine dorsal/ventral from filename: %s\n', swrFile);
            continue;
        end

        % Load SWR intervals
        load(swrFile);
        if strcmp(swr_location, 'Dorsal')
            probeSWR = dorsal_swr;
        else
            probeSWR = ventral_swr;
        end
        probeSWR.probeID = probeID;
        probeSWR.iv      = ResizeIV(cfg_iv, probeSWR.iv);

        if length(probeSWR.iv.tstart) < SWR_MIN_COUNT
            fprintf('  [skip] Too few SWRs (%d) in %s\n', ...
                length(probeSWR.iv.tstart), swrFile);
            continue;
        end

        % ---- Find best SWR channel depth (CA1 band filter) ------------------
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

        % ---- Subset to CA1 band cells on this probe -------------------------
        probeIdx    = find(strcmp(S.usr.probeID, probeID));
        depthsProbe = S.usr.recording_depths(probeIdx);
        ca1Mask     = abs(depthsProbe - swrChDepth) <= CA1_DEPTH_THRESH_UM;
        ca1Idx      = probeIdx(ca1Mask);

        if isempty(ca1Idx)
            fprintf('  [skip] No cells within CA1 depth band for %s\n', swrFile);
            continue;
        end

        fprintf('  [info] %d / %d probe cells pass CA1 depth filter\n', ...
            sum(ca1Mask), length(probeIdx));

        % ---- Decorrelation mask ---------------------------------------------
        decorr_keep = getDecorrMask(S, ca1Idx, probeID, ExpKeys);
        n_before    = numel(ca1Idx);
        ca1Idx      = ca1Idx(decorr_keep);
        n_removed   = n_before - numel(ca1Idx);
        if n_removed > 0
            fprintf('  [info] Decorr mask removed %d / %d CA1-band cells on %s\n', ...
                n_removed, n_before, probeID);
        end

        if isempty(ca1Idx)
            fprintf('  [skip] No cells remain after decorr mask for %s\n', swrFile);
            continue;
        end

        S_probe = subsetSpikes(S, ca1Idx);

        % Compute SWR participation with outcome-group labelling
        res = processCells(S_probe, probeSWR, sessionID, swr_location, ...
            ExpKeys.sessend, MIN_BASELINE_FR, odor_to_outcome);

        % Append to accumulators
        all_session           = [all_session;           res.session];
        all_label             = [all_label;             res.label];
        all_location          = [all_location;          res.location];
        all_cell_type         = [all_cell_type;         res.cell_type];
        all_outcome_group     = [all_outcome_group;     res.outcome_group];
        all_trial_type_set    = [all_trial_type_set;    res.trial_type_set];
        all_odor_outcome      = [all_odor_outcome;      res.odor_outcome];
        all_reward_dir        = [all_reward_dir;        res.reward_dir];
        all_airpuff_dir       = [all_airpuff_dir;       res.airpuff_dir];
        all_stimulus_encoding = [all_stimulus_encoding; res.stimulus_encoding];
        all_pct_participation = [all_pct_participation; res.pct_participation];
        all_swr_count         = [all_swr_count;         res.swr_count];

    end % SWR file loop
end % folder loop

%% =========================================================================
%% Build and save output table
%% =========================================================================

celltable = table( ...
    all_session,           ...
    all_label,             ...
    all_location,          ...
    all_cell_type,         ...
    all_outcome_group,     ...
    all_trial_type_set,    ...
    all_odor_outcome,      ...
    all_reward_dir,        ...
    all_airpuff_dir,       ...
    all_stimulus_encoding, ...
    all_pct_participation, ...
    all_swr_count,         ...
    'VariableNames', { ...
    'Session', 'Label', 'Location', 'CellType', ...
    'OutcomeGroup', 'TrialTypeSet', 'OdorOutcome', ...
    'RewardDirection', 'AirpuffDirection', ...
    'StimulusEncoding', ...
    'PctSWRParticipation', 'SWRCount'});

if isscalar(folder_list)
    filename     = char(strcat(ExpKeys.subject, '_', ExpKeys.date, ...
        '_SWR_participation_by_outcome.mat'));
    filename_csv = char(strcat(ExpKeys.subject, '_', ExpKeys.date, ...
        '_SWR_participation_by_outcome.csv'));
else
    folderPath   = 'D:\vCA1';
    baseName     = 'allSessions_SWR_participation_by_outcome';
    filename     = fullfile(folderPath, strcat(baseName, '.mat'));
    filename_csv = fullfile(folderPath, strcat(baseName, '.csv'));
end

save(filename, 'celltable', 'all_ExpKeys', 'cfg_iv');
writetable(celltable, filename_csv);
fprintf('\nResults saved to:\n  %s\n  %s\n', filename, filename_csv);

%% =========================================================================
%% Statistics
%% =========================================================================
% Three analyses, each with a distinct question:
%
%   Analysis 1 – WHICH OUTCOME drives differential SWR participation?
%     Restricted to stimulus-encoding neurons only (StimulusEncoding == 1).
%     Kruskal-Wallis across Reward / Airpuff / Neutral (3 groups, at most
%     3 post-hoc pairs). Unclassified neurons are excluded so the omnibus
%     test and Bonferroni correction are not inflated by an uninformative
%     catch-all group.
%
%   Analysis 2 – REWARD vs AIRPUFF (stimulus-encoding neurons only)
%     Direct pairwise Wilcoxon rank-sum. Membership defined by trial-type
%     SET containing the relevant outcome (broader than primary group), so
%     a neuron selective for a reward odor that also has an airpuff response
%     appears on both sides. Any such overlap is counted and flagged.
%
%   Analysis 3 – ENCODING vs UNCLASSIFIED
%     Wilcoxon rank-sum on the full population (all neurons passing the
%     depth / firing-rate filter). Reprises the comparison from v4 so the
%     broader context is preserved alongside the finer-grained analyses.

enc_mask  = all_stimulus_encoding == 1;   % stimulus-encoding neurons
unc_mask  = all_stimulus_encoding == 0;   % unclassified neurons

% Outcome groups for Analyses 1 & 2 (encoding neurons only)
enc_outcome_groups = {'Airpuff', 'Reward', 'Neutral'};
regions            = {'All', 'Dorsal', 'Ventral'};

% Helper: check whether trial-type set string contains a given outcome
% Uses exact token matching after splitting on commas (safer than strfind)
tt_contains = @(ttset_col, outcome) cellfun(@(s) any(strcmp(strsplit(s, ','), outcome)), ttset_col);

%% =========================================================================
%% Analysis 1: Kruskal-Wallis across Reward / Airpuff / Neutral
%%             (stimulus-encoding neurons only)
%% =========================================================================

fprintf('\n====================================================================\n');
fprintf('ANALYSIS 1: SWR Participation by Outcome Group\n');
fprintf('  Stimulus-encoding neurons only  |  Kruskal-Wallis + post-hoc\n');
fprintf('====================================================================\n');

for iReg = 1:length(regions)
    reg = regions{iReg};

    switch reg
        case 'All';     rmask = true(size(all_pct_participation));
        case 'Dorsal';  rmask = strcmp(all_location, 'Dorsal');
        case 'Ventral'; rmask = strcmp(all_location, 'Ventral');
    end

    base_mask = enc_mask & rmask;

    fprintf('\n--- %s CA1 ---\n', reg);

    group_data  = {};
    group_sizes = zeros(1, 3);
    kw_data     = [];
    kw_groups   = {};

    for iG = 1:length(enc_outcome_groups)
        grp  = enc_outcome_groups{iG};
        mask = base_mask & strcmp(all_outcome_group, grp);
        pct  = all_pct_participation(mask);

        group_data{iG}  = pct;
        group_sizes(iG) = numel(pct);

        fprintf('  %-10s  n = %3d,  median = %5.2f%%,  mean = %5.2f%%\n', ...
            grp, numel(pct), ...
            ternary(~isempty(pct), median(pct), NaN), ...
            ternary(~isempty(pct), mean(pct),   NaN));

        kw_data   = [kw_data;   pct(:)];                          %#ok<AGROW>
        kw_groups = [kw_groups; repmat({grp}, numel(pct), 1)];    %#ok<AGROW>
    end

    non_empty = find(group_sizes > 0);

    if numel(non_empty) >= 2 && numel(kw_data) >= 4
        [p_kw, tbl_kw, stats_kw] = kruskalwallis(kw_data, kw_groups, 'off');
        fprintf('  Kruskal-Wallis: p = %.4f,  chi2 = %.3f,  df = %d\n', ...
            p_kw, tbl_kw{2,5}, tbl_kw{2,3});
    else
        fprintf('  Kruskal-Wallis: insufficient data for test\n');
        continue;
    end

    % ---- Post-hoc: multcompare (Dunn's test with Bonferroni correction) -----
    % multcompare uses the rank sums from the KW stats struct directly,
    % avoiding the manual rank-averaging bug in the previous implementation.
    % Bonferroni correction is applied over the number of pairs actually
    % tested (adaptive, reflecting the tests genuinely conducted).

    pairs      = nchoosek(non_empty, 2);
    n_pairs    = size(pairs, 1);
    alpha_bonf = 0.05 / max(n_pairs, 1);

    [mc_results, ~] = multcompare(stats_kw, 'CType', 'bonferroni', 'Display', 'off');
    % mc_results columns: [groupA, groupB, lower_CI, mean_diff, upper_CI, p_value]

    fprintf('  Post-hoc Tukey-Kramer on ranks (Bonferroni alpha = %.4f, %d actual pairs):\n', ...
        alpha_bonf, n_pairs);

    for iPair = 1:size(mc_results, 1)
        iA = mc_results(iPair, 1);
        iB = mc_results(iPair, 2);

        % Skip if either group is empty in this region
        if group_sizes(iA) < 2 || group_sizes(iB) < 2
            continue;
        end

        gA   = enc_outcome_groups{iA};
        gB   = enc_outcome_groups{iB};
        pA   = group_data{iA};
        pB   = group_data{iB};
        nA   = numel(pA);
        nB   = numel(pB);
        p_mc = mc_results(iPair, 6);

        % Rank-biserial effect size via Mann-Whitney U
        [~, ~, stats_mw] = ranksum(pA, pB);
        U_val = stats_mw.ranksum - nA*(nA+1)/2;
        r_rb  = 1 - (2 * U_val) / (nA * nB);

        sig_str = ternary(p_mc < alpha_bonf, '  *', '');
        fprintf('    %-10s vs %-10s : p = %.4f,  r_rb = %+.3f%s\n', ...
            gA, gB, p_mc, r_rb, sig_str);
    end
end

%% =========================================================================
%% Analysis 2: Reward-encoding vs Airpuff-encoding
%%             (stimulus-encoding neurons only)
%% =========================================================================
% Membership: trial-type SET contains 'Reward' (or 'Airpuff'), restricted
% to encoding neurons. This is broader than the primary OutcomeGroup
% assignment, so a neuron with both a reward-odor preference and an airpuff
% response appears on both sides and is flagged below.

fprintf('\n====================================================================\n');
fprintf('ANALYSIS 2: Reward-encoding vs Airpuff-encoding Neurons\n');
fprintf('  Stimulus-encoding neurons only  |  Wilcoxon rank-sum\n');
fprintf('  (Membership: trial-type SET contains the relevant outcome)\n');
fprintf('====================================================================\n');

reward_mask  = enc_mask & tt_contains(all_trial_type_set, 'Reward');
airpuff_mask = enc_mask & tt_contains(all_trial_type_set, 'Airpuff');
overlap_n    = sum(reward_mask & airpuff_mask);

if overlap_n > 0
    fprintf('\n  Note: %d neuron(s) encode BOTH Reward and Airpuff outcomes\n', overlap_n);
    fprintf('        and appear on both sides of each comparison below.\n');
end

for iReg = 1:length(regions)
    reg = regions{iReg};

    switch reg
        case 'All';     rmask = true(size(all_pct_participation));
        case 'Dorsal';  rmask = strcmp(all_location, 'Dorsal');
        case 'Ventral'; rmask = strcmp(all_location, 'Ventral');
    end

    pct_reward  = all_pct_participation(rmask & reward_mask);
    pct_airpuff = all_pct_participation(rmask & airpuff_mask);

    fprintf('\n--- %s CA1 ---\n', reg);
    fprintf('  Reward-encoding   n = %3d,  median = %5.2f%%,  mean = %5.2f%%\n', ...
        numel(pct_reward), ...
        ternary(~isempty(pct_reward),  median(pct_reward),  NaN), ...
        ternary(~isempty(pct_reward),  mean(pct_reward),    NaN));
    fprintf('  Airpuff-encoding  n = %3d,  median = %5.2f%%,  mean = %5.2f%%\n', ...
        numel(pct_airpuff), ...
        ternary(~isempty(pct_airpuff), median(pct_airpuff), NaN), ...
        ternary(~isempty(pct_airpuff), mean(pct_airpuff),   NaN));

    if numel(pct_reward) >= 2 && numel(pct_airpuff) >= 2
        [p_ra, ~, stats_ra] = ranksum(pct_reward, pct_airpuff);
        if isfield(stats_ra, 'zval')
            fprintf('  Wilcoxon rank-sum: p = %.4f,  z = %.3f\n', p_ra, stats_ra.zval);
        else
            fprintf('  Wilcoxon rank-sum: p = %.4f\n', p_ra);
        end
    else
        fprintf('  Wilcoxon rank-sum: insufficient data (n_reward=%d, n_airpuff=%d)\n', ...
            numel(pct_reward), numel(pct_airpuff));
    end
end

%% =========================================================================
%% Analysis 3: Stimulus-encoding vs Unclassified  (full population)
%% =========================================================================
% Wilcoxon rank-sum on all neurons passing depth / firing-rate filters.
% Preserves the contrast from swr_participation_stimulus_encoding_v4.m
% so the broader context is available alongside the finer-grained analyses.

fprintf('\n====================================================================\n');
fprintf('ANALYSIS 3: Stimulus-Encoding vs Unclassified Neurons\n');
fprintf('  Full population  |  Wilcoxon rank-sum\n');
fprintf('====================================================================\n');

for iReg = 1:length(regions)
    reg = regions{iReg};

    switch reg
        case 'All';     rmask = true(size(all_pct_participation));
        case 'Dorsal';  rmask = strcmp(all_location, 'Dorsal');
        case 'Ventral'; rmask = strcmp(all_location, 'Ventral');
    end

    pct_enc   = all_pct_participation(rmask & enc_mask);
    pct_unenc = all_pct_participation(rmask & unc_mask);

    fprintf('\n--- %s CA1 ---\n', reg);
    fprintf('  Encoding      n = %3d,  median = %5.2f%%,  mean = %5.2f%%\n', ...
        numel(pct_enc), ...
        ternary(~isempty(pct_enc),   median(pct_enc),   NaN), ...
        ternary(~isempty(pct_enc),   mean(pct_enc),     NaN));
    fprintf('  Unclassified  n = %3d,  median = %5.2f%%,  mean = %5.2f%%\n', ...
        numel(pct_unenc), ...
        ternary(~isempty(pct_unenc), median(pct_unenc), NaN), ...
        ternary(~isempty(pct_unenc), mean(pct_unenc),   NaN));

    if numel(pct_enc) >= 2 && numel(pct_unenc) >= 2
        [p_eu, ~, stats_eu] = ranksum(pct_enc, pct_unenc);
        if isfield(stats_eu, 'zval')
            fprintf('  Wilcoxon rank-sum: p = %.4f,  z = %.3f\n', p_eu, stats_eu.zval);
        else
            fprintf('  Wilcoxon rank-sum: p = %.4f\n', p_eu);
        end
    else
        fprintf('  Wilcoxon rank-sum: insufficient data\n');
    end
end

%% =========================================================================
%% Summary: neurons participating in ≥10% of SWRs, by outcome group
%% =========================================================================
% Encoding neurons only, matching the scope of Analyses 1 & 2.

PCT_THRESH = 10;

fprintf('\n====================================================================\n');
fprintf('Neurons participating in ≥%d%% of SWRs, by outcome group\n', PCT_THRESH);
fprintf('  (Stimulus-encoding neurons only)\n');
fprintf('====================================================================\n');

for iReg = 1:length(regions)
    reg = regions{iReg};
    switch reg
        case 'All';     rmask = true(size(all_pct_participation));
        case 'Dorsal';  rmask = strcmp(all_location, 'Dorsal');
        case 'Ventral'; rmask = strcmp(all_location, 'Ventral');
    end

    fprintf('\n--- %s CA1 ---\n', reg);
    for iG = 1:length(enc_outcome_groups)
        grp  = enc_outcome_groups{iG};
        mask = enc_mask & rmask & strcmp(all_outcome_group, grp);
        pct  = all_pct_participation(mask);
        if isempty(pct)
            fprintf('  %-10s : no neurons\n', grp);
            continue;
        end
        n_total  = numel(pct);
        n_active = sum(pct >= PCT_THRESH);
        fprintf('  %-10s : %d / %d (%.1f%%) participate in ≥%d%% of SWRs\n', ...
            grp, n_active, n_total, 100*n_active/n_total, PCT_THRESH);
    end
end

%% =========================================================================
%% Bar Plot: Mean SWR Participation by Outcome Group and Region
%% =========================================================================
% Shows mean +/- SEM participation rate for stimulus-encoding neurons
% grouped by primary outcome (Reward / Airpuff / Neutral), split by
% Dorsal vs Ventral CA1.
%
% Layout: 1 x 2 subplots (Dorsal | Ventral), three grouped bars each.
% Individual neuron data points are overlaid as a scatter to show
% distribution shape alongside the mean.

plot_groups  = {'Airpuff', 'Reward', 'Neutral'};
plot_regions = {'Dorsal', 'Ventral'};
nGroups      = length(plot_groups);
nRegions     = length(plot_regions);

% Colours: Reward=blue, Airpuff=red, Neutral=grey
group_colors = [ ...
    0.00  0.45  0.74;   % Airpuff  - muted blue
    0.85  0.33  0.10;   % Reward - muted orange
    0.93  0.69  0.13;   % Neutral - yellow
];

% Pre-compute means, SEMs, and raw data per region x group
means = nan(nRegions, nGroups);
sems  = nan(nRegions, nGroups);
raw   = cell(nRegions, nGroups);

for iReg = 1:nRegions
    reg   = plot_regions{iReg};
    rmask = strcmp(all_location, reg);

    for iG = 1:nGroups
        grp   = plot_groups{iG};
        % Encoding neurons with this primary outcome group only
        gmask = rmask & enc_mask & strcmp(all_outcome_group, grp);

        pct = all_pct_participation(gmask);
        raw{iReg, iG} = pct;

        if ~isempty(pct)
            means(iReg, iG) = mean(pct);
            sems(iReg, iG)  = std(pct) / sqrt(numel(pct));
        end
    end
end

% ---- Plot -------------------------------------------------------------------
fig_bar = figure('Name', 'SWR Participation by Outcome Group & Region', ...
    'Position', [100 100 860 500], 'Color', 'w');

bar_width  = 0.5;    % width of each individual bar
group_gap  = 1.0;     % spacing between outcome groups
jitter_amt = 0.12;    % horizontal jitter for scatter

for iReg = 1:nRegions
    ax = subplot(1, 2, iReg);
    hold(ax, 'on');

    % X positions: one cluster per group
    x_centers = (1:nGroups) * group_gap;

    for iG = 1:nGroups
        xc  = x_centers(iG);
        clr = group_colors(iG, :);
        m   = means(iReg, iG);
        s   = sems(iReg, iG);

        if isnan(m);  continue;  end

        % Bar
        bar(ax, xc, m, bar_width, ...
            'FaceColor', clr, 'EdgeColor', 'none', 'FaceAlpha', 0.88);

        % Error bar (SEM)
        errorbar(ax, xc, m, s, ...
            'Color', [0.25 0.25 0.25], ...   % CHANGED: was 'k'
            'LineWidth', 1.0, ...             % CHANGED: was 2
            'CapSize', 0, ...                 % CHANGED: was 8
            'LineStyle', 'none');

        % Scatter overlay (individual neurons)
        pct = raw{iReg, iG};
        if ~isempty(pct)
            jitter = (rand(numel(pct), 1) - 0.5) * jitter_amt * 2;
            scatter(ax, xc + jitter, pct, 30, ...
                [0.92 0.92 0.92], 'filled', ...      % CHANGED: light gray instead of clr
                'MarkerFaceAlpha', 0.85, ...          % CHANGED: more opaque
                'MarkerEdgeColor', [0.3 0.3 0.3], ...% CHANGED: add thin dark edge back
                'LineWidth', 0.6);                    % NEW: thin edge line
        end
    end

    % Formatting
    set(ax, 'XTick', x_centers, 'XTickLabel', plot_groups, ...
        'FontSize', 14, ...
        'FontName', 'Helvetica', ...          % NEW
        'TickDir', 'out', ...
        'Box', 'off', ...
        'LineWidth', 1.0, ...                 % NEW
        'XColor', [0.2 0.2 0.2], ...          % NEW
        'YColor', [0.2 0.2 0.2]);             % NEW
    xtickangle(ax, 0);
    xlim(ax, [x_centers(1) - 0.5, x_centers(end) + 0.5]);
    ylim_auto = ylim(ax);
    ylim(ax, [0, max(ylim_auto(2), 1) * 1.12]);   % CHANGED: added *1.12 headroom
    ylabel(ax, 'SWR Participation (%)');
    title(ax, sprintf('%s CA1', plot_regions{iReg}), ...
        'FontWeight', 'normal', ...   % CHANGED: was 'bold'
        'FontSize', 15, ...           % NEW
        'Color', [0.15 0.15 0.15]);   % NEW

    % Annotate n per bar below x-axis
    for iG = 1:nGroups
        n_neurons = numel(raw{iReg, iG});
        if n_neurons == 0;  continue;  end
        xc = x_centers(iG);
        yl = ylim(ax);
        t = text(ax, xc, yl(1) - diff(yl)*0.1, sprintf('n=%d', n_neurons), ...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'top', ...
            'FontSize', 10, 'Color', [0.5 0.5 0.5]);
        t.Clipping = 'off';
    end

    hold(ax, 'off');
end

% Save figure as vector PDF alongside the data file
if isscalar(folder_list)
    fig_out = char(strcat(ExpKeys.subject, '_', ExpKeys.date, ...
        '_SWR_participation_barplot.png'));
else
    fig_out = fullfile('D:\vCA1', 'allSessions_SWR_participation_barplot.png');
end
exportgraphics(fig_bar, fig_out, 'Resolution', 600);
fprintf('\nBar plot saved to: %s\n', fig_out);

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function v = ternary(cond, a, b)
% ternary  Inline if-else helper for scalar display values.
if cond;  v = a;  else;  v = b;  end
end

% -------------------------------------------------------------------------

function csc = loadLFP(cfg, cfg_md, ExpKeys)
% loadLFP  Load LFP for the requested probe, apply true-depth correction,
%          and sort channels dorsal (shallow) -> ventral (deep).

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

% -------------------------------------------------------------------------

function lfp = selectSWRLfp(csc, ExpKeys, cfg)
% selectSWRLfp  Extract the SWR channel for the chosen probe/region.
% cfg.probe (1 = imec0, 2 = imec1) and cfg.region ('dorsal'/'ventral')
% determine which ExpKeys field to read.

probeStr = sprintf('imec%d', cfg.probe - 1);   % 1->imec0, 2->imec1
regionStr_lower = lower(cfg.region);
fieldName = sprintf('%s_%s_SWR_channel', probeStr, regionStr_lower);

if ~isfield(ExpKeys, fieldName)
    % Graceful fallback: warn and use imec0 equivalent
    warning('selectSWRLfp: ExpKeys field ''%s'' not found; falling back to imec0_%s_SWR_channel.', ...
        fieldName, regionStr_lower);
    fieldName = sprintf('imec0_%s_SWR_channel', regionStr_lower);
end
swrCh = ExpKeys.(fieldName);

swrIdx             = find(cellfun(@(x) strcmp(x, swrCh), csc.label), 1);
lfp.SWR            = csc;
lfp.SWR.data       = lfp.SWR.data(swrIdx, :);
lfp.SWR.cfg.depths = lfp.SWR.cfg.depths(swrIdx);
lfp.SWR.label      = lfp.SWR.label(swrIdx);
end

% -------------------------------------------------------------------------

function keep = getDecorrMask(S, idx, probeID, ExpKeys)
% getDecorrMask  Return logical vector; TRUE = not on a decorrelated shank.

n    = numel(idx);
keep = true(n, 1);

for iCell = 1:n
    shankID   = S.usr.shank_ids(idx(iCell));
    fieldName = sprintf('%s_shank%d_decorr_times', probeID, shankID);

    if ~isfield(ExpKeys, fieldName);  continue;  end

    val = ExpKeys.(fieldName);
    if iscell(val);  val = val{1};  end
    val = double(val);

    if val >= 0
        keep(iCell) = false;
    end
end
end

% -------------------------------------------------------------------------

function S_sub = subsetSpikes(S, idx)
% subsetSpikes  Return a copy of S containing only the cells at indices idx.

S_sub       = S;
S_sub.label = S.label(idx);
S_sub.t     = S.t(idx);
f = fieldnames(S.usr);
for iF = 1:length(f)
    val = S.usr.(f{iF});
    if iscell(val) && numel(val) == numel(S.label)
        S_sub.usr.(f{iF}) = val(idx);
    elseif isnumeric(val) && numel(val) == numel(S.label)
        S_sub.usr.(f{iF}) = val(idx);
    end
end
end

% -------------------------------------------------------------------------

function res = processCells(S, probeSWR, sessionID, swr_location, ...
    sessend, min_baseline_fr, odor_to_outcome)
% processCells  Compute SWR participation % and outcome-group label for each
%               qualifying cell in swr_location.
%
% Trial-type logic mirrors computeCoactivity() in swr_coactivity_same_vs_diff_v2.m:
%   - Trial-type SET = {odor outcome (if defined)} ∪ {'Reward' if reward_dir ~= 'None'}
%                                                   ∪ {'Airpuff' if airpuff_dir ~= 'None'}
%   - StimulusEncoding = true iff at least one excitatory/increase response
%   - Primary OutcomeGroup assigned by priority: Reward > Airpuff > Neutral > Unclassified
%
% ALL cells passing the depth/firing-rate filter are included (not just
% stimulus-encoding ones) so Unclassified neurons appear as a comparison group.

nSWR = length(probeSWR.iv.tstart);

res.session           = {};
res.label             = {};
res.location          = {};
res.cell_type         = {};
res.outcome_group     = {};
res.trial_type_set    = {};
res.odor_outcome      = {};
res.reward_dir        = {};
res.airpuff_dir       = {};
res.stimulus_encoding = [];
res.pct_participation = [];
res.swr_count         = [];

for iCell = 1:length(S.label)

    % Region filter
    if ~strcmpi(S.usr.HC_location{iCell}, swr_location)
        continue;
    end

    % Firing-rate filter
    baseline_fr = numel(S.t{iCell}) / sessend;
    if baseline_fr < min_baseline_fr
        continue;
    end

    % ---- Retrieve selectivity fields ------------------------------------
    odor_val        = S.usr.odor{iCell};
    odor_dir_val    = S.usr.odor_dir{iCell};
    reward_dir_val  = S.usr.reward_dir{iCell};
    airpuff_dir_val = S.usr.airpuff_dir{iCell};

    % ---- Excitatory/increase filter (mirrors coactivity script) ---------
    has_odor_exc    = strcmpi(odor_dir_val,    'excitatory') || strcmpi(odor_dir_val,    'increase');
    has_reward_inc  = strcmpi(reward_dir_val,  'increase');
    has_airpuff_inc = strcmpi(airpuff_dir_val, 'increase');
    stimulus_encoding = has_odor_exc || has_reward_inc || has_airpuff_inc;

    % ---- Odor identity / outcome ----------------------------------------
    odor_is_defined = ~strcmpi(odor_val, 'None') && ~isempty(odor_val);

    if odor_is_defined && isKey(odor_to_outcome, odor_val)
        odor_outcome_val = odor_to_outcome(odor_val);
    else
        odor_outcome_val = 'None';
    end

    % ---- Trial-type set (same logic as computeCoactivity) ---------------
    % Include odor outcome if the odor is defined,
    % plus direct reward/airpuff responses. Given every neuron is at
    % least excitatory for one thing, this helps to push neurons to be
    % labelled as the thing they're most responsive for
    tt = {};
    if odor_is_defined && ~strcmp(odor_outcome_val, 'None') && ~strcmp(odor_outcome_val, 'inhibitory') ...
            && ~strcmp(odor_outcome_val, 'decrease')
        tt{end+1} = odor_outcome_val;
    end
    if ~strcmpi(reward_dir_val,  'None') && ~isempty(reward_dir_val) && ~strcmpi(reward_dir_val,  'decrease')
        tt{end+1} = 'Reward';
    end
    if ~strcmpi(airpuff_dir_val, 'None') && ~isempty(airpuff_dir_val) && ~strcmpi(airpuff_dir_val,  'decrease')
        tt{end+1} = 'Airpuff';
    end
    tt = unique(tt);

    % ---- Primary outcome group (priority: Reward > Neutral > Airpuff) ---
    if any(strcmp(tt, 'Reward'))
        outcome_group = 'Reward';
    elseif any(strcmp(tt, 'Neutral'))
        outcome_group = 'Neutral';
    elseif any(strcmp(tt, 'Airpuff'))
        outcome_group = 'Airpuff';
    else
        outcome_group = 'Unclassified';
    end

    % Store trial-type set as comma-separated string for table readability
    if isempty(tt)
        tt_str = 'None';
    else
        tt_str = strjoin(sort(tt), ',');
    end

    % ---- Cell-type label (for downstream compatibility) -----------------
    is_odor_selective     = odor_is_defined;
    is_outcome_responsive = (~strcmpi(reward_dir_val,  'None') && ~isempty(reward_dir_val)) || ...
        (~strcmpi(airpuff_dir_val, 'None') && ~isempty(airpuff_dir_val));
    if is_odor_selective && is_outcome_responsive
        cell_type = 'Odor_and_Outcome';
    elseif is_odor_selective
        cell_type = 'Odor_Selective';
    elseif is_outcome_responsive
        cell_type = 'Outcome_Responsive';
    else
        cell_type = 'Unclassified';
    end

    % ---- Compute % SWR participation ------------------------------------
    currCell_t               = S.t{iCell};
    swr_participation_binary = zeros(1, nSWR);

    for iSWR = 1:nSWR
        t_start = probeSWR.iv.tstart(iSWR);
        t_end   = probeSWR.iv.tend(iSWR);
        if sum(currCell_t > t_start & currCell_t < t_end) > 0
            swr_participation_binary(iSWR) = 1;
        end
    end

    pct_participation = mean(swr_participation_binary) * 100;

    % ---- Append to result struct ----------------------------------------
    res.session           = [res.session;           {sessionID}];
    res.label             = [res.label;             {S.label{iCell}}];
    res.location          = [res.location;          {S.usr.HC_location{iCell}}];
    res.cell_type         = [res.cell_type;         {cell_type}];
    res.outcome_group     = [res.outcome_group;     {outcome_group}];
    res.trial_type_set    = [res.trial_type_set;    {tt_str}];
    res.odor_outcome      = [res.odor_outcome;      {odor_outcome_val}];
    res.reward_dir        = [res.reward_dir;        {reward_dir_val}];
    res.airpuff_dir       = [res.airpuff_dir;       {airpuff_dir_val}];
    res.stimulus_encoding = [res.stimulus_encoding; stimulus_encoding];
    res.pct_participation = [res.pct_participation; pct_participation];
    res.swr_count         = [res.swr_count;         nSWR];

end % cell loop
end % processCells