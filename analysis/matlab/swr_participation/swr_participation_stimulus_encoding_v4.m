%% Sharp Wave Ripple Participation - Stimulus-Encoding vs. Unclassified
% Tests whether stimulus-encoding neurons participate in a higher percentage
% of SWRs than unclassified neurons, using a Wilcoxon rank-sum test.
% Analysis is run overall and separately for dorsal and ventral CA1.
%
% Grouping:
%   StimulusEncoding = true  : has at least one EXCITATORY / INCREASE response
%                              (odor_dir == 'excitatory' or 'increase', OR
%                               reward_dir == 'increase', OR
%                               airpuff_dir == 'increase')
%   StimulusEncoding = false : Unclassified (no excitatory/increase response)
%
% Key design decisions:
%   - Each neuron is evaluated only against SWRs from its own region
%     (dorsal cells -> dorsal SWRs, ventral cells -> ventral SWRs).
%   - CA1 depth band filter (±CA1_DEPTH_THRESH_UM around best SWR channel)
%     excludes cells too far from the CA1 recording site.
%   - Participation = fired at least one spike during a SWR window.
%   - Comparison metric: % of SWRs a neuron participates in.
%   - Statistical test: Wilcoxon rank-sum (ranksum), two-sided.
%
% EG May 2026 - stripped from swr_participation_stimulus_encoding_v1.m
%
% DEPENDENCIES
%   mvdmlab / odor-pixels codebase functions:
%     LoadExpKeys, FindFiles, load_NPX_spike, ResizeIV, restrict, iv,
%     SelectTS, selectSWRLfp (or equivalent)
%   Data files per session folder:
%     clean_units_imec0.mat                    (required)
%     clean_units_imec1.mat                    (optional, dual-probe)
%     selective_permissive_imec0_all.mat       (optional; absent -> Unclassified)
%     selective_permissive_imec1_all.mat       (optional)
%     imec0_clean_lfp.mat                      (required, for SWR channel depth)
%     imec1_clean_lfp.mat                      (optional)
%     *SWR.mat                                 (probe ID + 'dorsal'/'ventral' in name)
%     ExpKeys file

%% -------------------------------------------------------------------------
%% Parameters
%% -------------------------------------------------------------------------

folder_list = {'D:\vCA1\M656\preprocessed\M656_2026_03_25'
               'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
               'D:\vCA1\M649\preprocessed\M649_2026_04_01'
               'D:\vCA1\M593\preprocessed\M593_2025_12_26'
               'D:\vCA1\M593\preprocessed\M593_2025_12_28'
               'D:\vCA1\M653\preprocessed\M653_2026_02_04'
               };

SWR_MIN_COUNT       = 20;    % min SWRs required to include a session/file
CA1_DEPTH_THRESH_UM = 300;   % ±µm around best SWR channel to define CA1 band

MIN_BASELINE_FR = 0;   % Hz; neurons below this across the whole session are excluded - consistent with encoding neurons 
EXCLUDE_NO_INCREASE = true; % set to true if you want only neurons that increase their firing rate to at least one stimulus to be used

spikeFileBase = 'selective_permissive';

cfg_iv   = [];
cfg_iv.d = [0, 0];   % expand each SWR interval (s) on each side

%% -------------------------------------------------------------------------
%% Initialize accumulators
%% -------------------------------------------------------------------------

all_session           = {};
all_label             = {};
all_location          = {};
all_cell_type         = {};
all_stimulus_encoding = [];   % logical: has ≥1 excitatory/increase response
all_pct_participation = [];
all_swr_count         = [];
all_ExpKeys           = {};

%% -------------------------------------------------------------------------
%% Main loop - folders
%% -------------------------------------------------------------------------

for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});

    if strcmp(folder_list{iFolder,1}, 'D:\vCA1\M656\preprocessed\M656_2026_03_25')
        CA1_DEPTH_THRESH_UM = 400;
        fprintf('increased CA1 depth threshold for M656_2026_03_25')
    else
        CA1_DEPTH_THRESH_UM = 300;
    end

    rng(2024);
    LoadExpKeys;
    all_ExpKeys{end+1} = ExpKeys;

    % ---- Depth calibration (probe 1) ----------------------------------------
    record_depth    = ExpKeys.probe1_Depth;
    cfg_md.maxdepth = 1000 * record_depth;

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

    % ---- Initialise selectivity fields to 'None' for every cell -------------
    nCells = length(S.label);
    S.usr.odor        = repmat({'None'}, 1, nCells);
    S.usr.odor_dir    = repmat({'None'}, 1, nCells);
    S.usr.reward_dir  = repmat({'None'}, 1, nCells);
    S.usr.airpuff_dir = repmat({'None'}, 1, nCells);

    % Build a label -> index map for fast lookup
    allUnitsMap = containers.Map(S.label, num2cell(1:nCells));

    % Merge selectivity from selective_permissive files (one per probe)
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

    % ---- Dorsal / ventral classification (per-probe threshold) --------------
    DH_threshold_imec0 = 3000;
    DH_threshold_imec1 = 300;

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

        % ---- Find best SWR channel depth (used for CA1 band filter only) ----
        cfg_lfp            = [];
        cfg_lfp.probe      = probeNum;
        cfg_lfp.region     = regionStr;

        csc        = loadLFP(cfg_lfp, cfg_md, ExpKeys);
        csc        = restrict(csc, 0, ExpKeys.sessend);
        lfp        = selectSWRLfp(csc, ExpKeys, cfg_lfp);
        swrChDepth = lfp.SWR.cfg.depths;   % scalar: true depth (µm) of best SWR channel

        fprintf('  [info] %s | best SWR ch depth = %.1f µm  (CA1 band: %.1f – %.1f µm)\n', ...
                swrFile, swrChDepth, ...
                swrChDepth - CA1_DEPTH_THRESH_UM, ...
                swrChDepth + CA1_DEPTH_THRESH_UM);

        % ---- Subset to cells on this probe within the CA1 band --------------
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

        % ---- Decorrelation mask: exclude cells on shanks with a decorr event --
        % Shanks with a non-negative decorr time in ExpKeys are considered
        % corrupted after that point and must be excluded entirely.
        decorr_keep = getDecorrMask(S, ca1Idx, probeID, ExpKeys);
        n_before     = numel(ca1Idx);
        ca1Idx       = ca1Idx(decorr_keep);
        n_removed    = n_before - numel(ca1Idx);
        if n_removed > 0
            fprintf('  [info] Decorr mask removed %d / %d CA1-band cells on %s\n', ...
                    n_removed, n_before, probeID);
        end

        if isempty(ca1Idx)
            fprintf('  [skip] No cells remain after decorr mask for %s\n', swrFile);
            continue;
        end

        S_probe = subsetSpikes(S, ca1Idx);

        % Compute SWR participation for each qualifying cell
        res = processCells(S_probe, probeSWR, sessionID, swr_location, ...
            ExpKeys.sessend, MIN_BASELINE_FR, EXCLUDE_NO_INCREASE);

        % Append to accumulators
        all_session           = [all_session;           res.session];
        all_label             = [all_label;             res.label];
        all_location          = [all_location;          res.location];
        all_cell_type         = [all_cell_type;         res.cell_type];
        all_stimulus_encoding = [all_stimulus_encoding; res.stimulus_encoding];
        all_pct_participation = [all_pct_participation; res.pct_participation];
        all_swr_count         = [all_swr_count;         res.swr_count];

    end % SWR file loop
end % folder loop

%% -------------------------------------------------------------------------
%% Build and save output table
%% -------------------------------------------------------------------------

celltable = table( ...
    all_session,           ...
    all_label,             ...
    all_location,          ...
    all_cell_type,         ...
    all_stimulus_encoding, ...
    all_pct_participation, ...
    all_swr_count,         ...
    'VariableNames', { ...
        'Session', 'Label', 'Location', 'CellType', ...
        'StimulusEncoding', ...
        'PctSWRParticipation', 'SWRCount'});

if isscalar(folder_list)
    filename     = char(strcat(ExpKeys.subject, '_', ExpKeys.date, ...
                               '_SWR_participation_stimenc.mat'));
    filename_csv = char(strcat(ExpKeys.subject, '_', ExpKeys.date, ...
                               '_SWR_participation_stimenc.csv'));
else
    folderPath   = 'D:\vCA1';
    baseName     = 'allSessions_SWR_participation_stimenc_inc';
    filename     = fullfile(folderPath, strcat(baseName, '.mat'));
    filename_csv = fullfile(folderPath, strcat(baseName, '.csv'));
end

save(filename, 'celltable', 'all_ExpKeys', 'cfg_iv');
writetable(celltable, filename_csv);

fprintf('\nResults saved to:\n  %s\n  %s\n', filename, filename_csv);

%% =========================================================================
%% Statistics: Wilcoxon rank-sum, stimulus-encoding vs. unclassified
%% =========================================================================

enc_mask   = all_stimulus_encoding == 1;
unenc_mask = all_stimulus_encoding == 0;

pct_enc   = all_pct_participation(enc_mask);
pct_unenc = all_pct_participation(unenc_mask);

fprintf('\n========== SWR Participation: Stimulus-Encoding vs. Unclassified ===========\n');
fprintf('  All regions combined\n');
fprintf('    Encoding   : n = %d,  median = %.2f%%,  mean = %.2f%%\n', numel(pct_enc),   median(pct_enc),   mean(pct_enc));
fprintf('    Unclassified: n = %d,  median = %.2f%%,  mean = %.2f%%\n', numel(pct_unenc), median(pct_unenc), mean(pct_unenc));
[p_all, h_all, stats_all] = ranksum(pct_enc, pct_unenc);
fprintf('    Wilcoxon rank-sum: p = %.4f,  z = %.3f\n', p_all, stats_all.zval);

% ---- Dorsal -----------------------------------------------------------------
d_enc   = all_pct_participation(enc_mask   & strcmp(all_location, 'Dorsal'));
d_unenc = all_pct_participation(unenc_mask & strcmp(all_location, 'Dorsal'));

fprintf('\n  Dorsal CA1\n');
fprintf('    Encoding   : n = %d,  median = %.2f%%, mean = %.2f%%\n',  numel(d_enc),   median(d_enc), mean(d_enc));
fprintf('    Unclassified: n = %d,  median = %.2f%%, mean = %.2f%%\n', numel(d_unenc), median(d_unenc), mean(d_unenc));
[p_d, h_d, stats_d] = ranksum(d_enc, d_unenc);
fprintf('    Wilcoxon rank-sum: p = %.4f,  z = %.3f\n', p_d, stats_d.zval);

% ---- Ventral ----------------------------------------------------------------
v_enc   = all_pct_participation(enc_mask   & strcmp(all_location, 'Ventral'));
v_unenc = all_pct_participation(unenc_mask & strcmp(all_location, 'Ventral'));

fprintf('\n  Ventral CA1\n');
fprintf('    Encoding   : n = %d,  median = %.2f%%, mean = %.2f%%\n', numel(v_enc),   median(v_enc), mean(v_enc));
fprintf('    Unclassified: n = %d,  median = %.2f%%, mean = %.2f%%\n', numel(v_unenc), median(v_unenc), mean(v_unenc));
[p_v, h_v, stats_v] = ranksum(v_enc, v_unenc);
fprintf('    Wilcoxon rank-sum: p = %.4f,  z = %.3f\n', p_v, stats_v.zval);

%% =========================================================================
%% Neurons active in ≥10% of SWRs
%% =========================================================================
% For each group (encoding / unclassified) and region, report what
% fraction of neurons participate in at least 10% of their local SWRs.

PCT_THRESH = 10;   % % SWR participation threshold

fprintf('\n=== Neurons participating in ≥%d%% of SWRs ===\n', PCT_THRESH);

groups      = {'Encoding', 'Unclassified'};
group_masks = {enc_mask,    unenc_mask};
regions     = {'All', 'Dorsal', 'Ventral'};

for iGroup = 1:2
    gmask = group_masks{iGroup};
    fprintf('\n  %s neurons\n', groups{iGroup});

    for iReg = 1:length(regions)
        reg = regions{iReg};

        switch reg
            case 'All'
                rmask = true(size(gmask));
            case 'Dorsal'
                rmask = strcmp(all_location, 'Dorsal');
            case 'Ventral'
                rmask = strcmp(all_location, 'Ventral');
        end

        pct = all_pct_participation(gmask & rmask);
        if isempty(pct)
            fprintf('    %-8s: no neurons\n', reg);
            continue;
        end

        n_total  = numel(pct);
        n_active = sum(pct >= PCT_THRESH);
        frac_pct = 100 * n_active / n_total;

        fprintf('    %-8s: %d / %d neurons (%.1f%%) participate in ≥%d%% of SWRs\n', ...
                reg, n_active, n_total, frac_pct, PCT_THRESH);
    end
end

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

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

% -------------------------------------------------------------------------

function lfp = selectSWRLfp(csc, ExpKeys, cfg)
% selectSWRLfp  Extract the SWR channel for the chosen probe region.
%
%   cfg.region     : 'dorsal' or 'ventral'

    if strcmpi(cfg.region, 'dorsal')
        swrCh  = ExpKeys.imec0_dorsal_SWR_channel;
    else
        swrCh  = ExpKeys.imec0_ventral_SWR_channel;
    end

    swrIdx             = find(cellfun(@(x) strcmp(x, swrCh),  csc.label), 1);
    lfp.SWR            = csc;
    lfp.SWR.data       = lfp.SWR.data(swrIdx, :);
    lfp.SWR.cfg.depths = lfp.SWR.cfg.depths(swrIdx);
    lfp.SWR.label      = lfp.SWR.label(swrIdx);
end

% -------------------------------------------------------------------------

function keep = getDecorrMask(S, idx, probeID, ExpKeys)
% getDecorrMask  Return a logical vector (length = numel(idx)) that is TRUE
%               for cells that are NOT on a decorrelated shank.
%
% A shank is considered decorrelated if ExpKeys contains a field named
%   <probeID>_shank<N>_decorr_times
% whose value is a cell scalar (or numeric scalar) with a non-negative number.
% Cells whose S.usr.shank_ids entry matches such a shank are excluded.
%
% Inputs:
%   S       - full spike struct (before subsetting)
%   idx     - indices into S of the cells to evaluate
%   probeID - 'imec0' or 'imec1'
%   ExpKeys - session ExpKeys struct

    n    = numel(idx);
    keep = true(n, 1);

    for iCell = 1:n
        shankID = S.usr.shank_ids(idx(iCell));   % scalar integer (0-based shank index)
        fieldName = sprintf('%s_shank%d_decorr_times', probeID, shankID);

        if ~isfield(ExpKeys, fieldName)
            continue;   % no decorr info for this shank -> keep
        end

        val = ExpKeys.(fieldName);
        if iscell(val)
            val = val{1};   % unwrap cell scalar
        end
        val = double(val);

        if val >= 0
            keep(iCell) = false;   % shank has a valid decorr event -> exclude
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

function res = processCells(S, probeSWR, sessionID, swr_location, sessend,...
    min_baseline_fr, exclude_no_increase)
% processCells  Compute SWR participation % for cells matching swr_location.
%
% Participation = fired at least one spike during a SWR window.
% Only cells whose HC_location matches swr_location are evaluated.
% Cells with a whole-session baseline firing rate < min_baseline_fr (Hz)
% are now fully excluded before any praticipation computation
%
% Inputs:
%   S            - spike struct (already subsetted to CA1 band)
%   probeSWR     - SWR struct with .iv.tstart / .iv.tend
%   sessionID    - string, e.g. 'M649-2026_03_30'
%   swr_location - 'Dorsal' or 'Ventral'
%   sess_end - scalar (s); end time of the recording session
%   min_baseline_fr - scalar (Hz); neurons firing below this are excluded 
%
% Output:
%   res - struct of column vectors / cell arrays, one row per qualifying cell

    nSWR = length(probeSWR.iv.tstart);

    res.session           = {};
    res.label             = {};
    res.location          = {};
    res.cell_type         = {};
    res.stimulus_encoding = [];
    res.pct_participation = [];
    res.swr_count         = [];

    for iCell = 1:length(S.label)

        % Only evaluate cells in the matching hippocampal region
        if ~strcmpi(S.usr.HC_location{iCell}, swr_location)
            continue;
        end

        % Baseline firing rate mask
        baseline_fr = numel(S.t{iCell}) / sessend;   % spikes / s over whole session
        if baseline_fr < min_baseline_fr
            continue;
        end

        % ---- Determine stimulus encoding status -----------------------------
        odor_val        = S.usr.odor{iCell};
        odor_dir_val    = S.usr.odor_dir{iCell};
        reward_dir_val  = S.usr.reward_dir{iCell};
        airpuff_dir_val = S.usr.airpuff_dir{iCell};

        % Toggle based on whether you want to include all stimulus
        % responsive neurons, or just the ones that increase their firing
        % rate to at least one stimulus 
        if exclude_no_increase
            % Excitatory / increase filter (mirrors swr_coactivity_same_vs_diff)
            has_odor_exc    = strcmpi(odor_dir_val,    'excitatory') || strcmpi(odor_dir_val,    'increase');
            has_reward_inc  = strcmpi(reward_dir_val,  'increase');
            has_airpuff_inc = strcmpi(airpuff_dir_val, 'increase');
    
            stimulus_encoding = has_odor_exc || has_reward_inc || has_airpuff_inc;
    
            % Convenience flags for cell_type label (direction-agnostic presence)
            is_odor_selective     = ~strcmpi(odor_val,        'None') && ~isempty(odor_val);
            is_outcome_responsive = (~strcmpi(reward_dir_val,  'None') && ~isempty(reward_dir_val)) || ...
                                    (~strcmpi(airpuff_dir_val, 'None') && ~isempty(airpuff_dir_val));
        else
            is_odor_selective     = ~strcmpi(odor_val,        'None') && ~isempty(odor_val);
            is_reward_responsive  = ~strcmpi(reward_dir_val,  'None') && ~isempty(reward_dir_val);
            is_airpuff_responsive = ~strcmpi(airpuff_dir_val, 'None') && ~isempty(airpuff_dir_val);
            is_outcome_responsive = is_reward_responsive || is_airpuff_responsive;

            stimulus_encoding = is_odor_selective || is_outcome_responsive;
        end


        % Four-way label retained for table inspection
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
        res.stimulus_encoding = [res.stimulus_encoding; stimulus_encoding];
        res.pct_participation = [res.pct_participation; pct_participation];
        res.swr_count         = [res.swr_count;         nSWR];

    end % cell loop
end % processCells