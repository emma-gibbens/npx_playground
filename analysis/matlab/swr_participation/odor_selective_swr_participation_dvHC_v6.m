%% Sharp Wave Ripple Participation Script - Dorsal vs Ventral Hippocampus
% Compares SWR participation by hippocampal location for ALL neurons,
% including odor-selective, outcome-responsive, and unclassified cells.
% Unclassified cells are loaded via load_NPX_spike('clean_units_imecX.mat').
% Selective cells are loaded from selective_imecX_all.mat (S_out) and their
% usr selectivity fields are merged into the all-units struct by label matching.
% EG May 2026 - adapted from odor_selective_swr_participation_dvHC_v4.m
%
% DEPENDENCIES
% - Requires multiple functions from the mvdmlab and odor-pixels codebase
% - Data files:
%   * clean_units_imec0.mat (required)
%   * clean_units_imec1.mat (optional, for dual-probe recordings)
%   * ExpKeys file (loaded via LoadExpKeys)
%   * selective_imec0_all.mat or selective_imec1_all.mat
%   * imec0_clean_lfp.mat (required)
%   * imec1_clean_lfp.mat (optional, for dual-probe recordings)

%% -------------------------------------------------------------------------
%% Parameters
%% -------------------------------------------------------------------------

folder_list = {'D:\vCA1\M656\preprocessed\M656_2026_03_25'};

SWR_MIN_COUNT          = 20;   % Minimum number of SWRs required to include a session
SWR_PARTICIPATION_THRESH = 10; % % of SWRs a cell must fire in to count as "participating"
CA1_DEPTH_THRESH_UM    = 400;  % Window (µm) around best SWR channel depth
                               % to define the CA1 band; cells outside this
                               % range are excluded from analysis.

cfg_iv   = [];
cfg_iv.d = [0, 0];    % Window (s) to expand each SWR interval by on each side

cfg_md   = [];                 % Will be populated per-folder below

%% -------------------------------------------------------------------------
%% Initialize accumulators
%% -------------------------------------------------------------------------

all_session           = {};
all_label             = {};
all_location          = {};
all_baseline_fr       = [];
all_mean_fr_allSWR    = [];
all_mean_fr_partSWR   = [];
all_pct_participation = [];
all_participates      = [];
all_swr_count         = [];
all_odor              = {};
all_odor_dir          = {};
all_odor_outcome      = {};
all_reward_dir        = {};
all_airpuff_dir       = {};
all_cell_type         = {};

%% -------------------------------------------------------------------------
%% Main loop - folders
%% -------------------------------------------------------------------------

for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});

    rng(2024);
    LoadExpKeys;

    % ---- Depth calibration (probe 1) ----------------------------------------
    record_depth        = ExpKeys.probe1_Depth;
    insertion_roll      = ExpKeys.probe1_insertion_roll;
    cfg_md.maxdepth     = 1000 * record_depth;

    % ---- Odor -> outcome map -------------------------------------------------
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

    % ---- Session duration (for baseline FR) ---------------------------------
    session_duration = ExpKeys.sessend - ExpKeys.sessstart;

    % ---- Load all clean units -----------------------------------------------
    % recording_depths is converted to true depth here at load time.
    % We do NOT redo the depth calculation later.
    S = load_NPX_spike('clean_units_imec0.mat');
    S.usr.recording_depths = (cfg_md.maxdepth - S.usr.recording_depths) ...
                              * cos(deg2rad(ExpKeys.probe1_insertion_roll));
    % Tag every imec0 cell so we can match it to imec0 SWRs later
    S.usr.probeID = repmat({'imec0'}, 1, length(S.label));

    if isfile('clean_units_imec1.mat')
        record_depth2 = ExpKeys.probe2_Depth;
        S2 = load_NPX_spike('clean_units_imec1.mat');
        S2.usr.recording_depths = (1000 * record_depth2 - S2.usr.recording_depths) ...
                                   * cos(deg2rad(ExpKeys.probe2_insertion_roll));
        % Tag imec1 cells
        S2.usr.probeID = repmat({'imec1'}, 1, length(S2.label));

        % Merge into S
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

    % ---- Merge selectivity fields from S_out structs ------------------------
    % Start every cell with 'None' for all selectivity fields.
    nCells = length(S.label);
    S.usr.odor        = repmat({'None'}, 1, nCells);
    S.usr.odor_dir    = repmat({'None'}, 1, nCells);
    S.usr.reward_dir  = repmat({'None'}, 1, nCells);
    S.usr.airpuff_dir = repmat({'None'}, 1, nCells);

    % Build a label->index map for the merged all-units struct
    allUnitsMap = containers.Map(S.label, num2cell(1:nCells));

    % For each selective file that exists, copy its usr fields into S
    selectiveFiles = {'selective_imec0_all.mat', 'selective_imec1_all.mat'};
    for iSF = 1:length(selectiveFiles)
        if exist(selectiveFiles{iSF}, 'file') ~= 2
            fprintf('  [info] %s not found - cells from this probe treated as Unclassified\n', ...
                    selectiveFiles{iSF});
            continue;
        end
        load(selectiveFiles{iSF});  % loads S_out
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
    % recording_depths is already true depth (converted at load time above).
    % DH threshold is per-probe. Used to label HC_location; CA1 band filtering
    % is done per SWR file below using the best SWR channel depth.
    DH_threshold_imec0 = 3000;
    DH_threshold_imec1 = 2500;

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

    % ---- SWR loop -----------------------------------------------------------
    % Each SWR file is probe-specific; we subset S to matching-probe cells
    % before calling processCells.
    AllSWR = FindFiles('*SWR.mat');
    sessionID = strcat(ExpKeys.subject, '-', ExpKeys.date);

    for currFile = 1:length(AllSWR)

        swrFile = AllSWR{currFile};

        % Identify which probe this SWR file belongs to
        if contains(swrFile, 'imec0')
            probeID  = 'imec0';
            probeNum = 1;
        elseif contains(swrFile, 'imec1')
            probeID  = 'imec1';
            probeNum = 2;
        else
            fprintf('  [skip] Cannot determine probe from filename: %s\n', swrFile);
            continue;
        end

        % Determine dorsal vs ventral region for this SWR file
        if contains(swrFile, 'dorsal')
            regionStr = 'dorsal';
        elseif contains(swrFile, 'ventral')
            regionStr = 'ventral';
        else
            fprintf('  [skip] Cannot determine dorsal/ventral from filename: %s\n', swrFile);
            continue;
        end

        % Load SWR intervals
        load(swrFile);
        if strcmp(regionStr, 'dorsal')
            probeSWR = dorsal_swr;
        else
            probeSWR = ventral_swr;
        end
        probeSWR.probeID = probeID;
        probeSWR.iv      = ResizeIV(cfg_iv, probeSWR.iv);

        % Skip sessions with too few SWRs
        if length(probeSWR.iv.tstart) < SWR_MIN_COUNT
            fprintf('  [skip] Too few SWRs (%d) in %s\n', length(probeSWR.iv.tstart), swrFile);
            continue;
        end

        % ---- Load LFP for this probe and find best SWR channel depth --------
        cfg_lfp        = [];
        cfg_lfp.probe  = probeNum;
        cfg_lfp.region = regionStr;
        cfg_lfp.useControl = false;

        csc = loadLFP(cfg_lfp, cfg_md, ExpKeys);
        csc = restrict(csc, 0, ExpKeys.sessend);

        lfp     = selectSWRLfp(csc, ExpKeys, cfg_lfp);
        swrChDepth = lfp.SWR.cfg.depths;   % scalar: true depth (µm) of best SWR channel

        fprintf('  [info] %s | best SWR ch depth = %.1f µm  (CA1 band: %.1f – %.1f µm)\n', ...
                swrFile, swrChDepth, ...
                swrChDepth - CA1_DEPTH_THRESH_UM, ...
                swrChDepth + CA1_DEPTH_THRESH_UM);

        % ---- Subset to probe cells that are within CA1 depth band -----------
        probeIdx  = find(strcmp(S.usr.probeID, probeID));
        depthsProbe = S.usr.recording_depths(probeIdx);
        ca1Mask   = abs(depthsProbe - swrChDepth) <= CA1_DEPTH_THRESH_UM;
        ca1Idx    = probeIdx(ca1Mask);

        if isempty(ca1Idx)
            fprintf('  [skip] No cells within CA1 depth band for %s\n', swrFile);
            continue;
        end

        fprintf('  [info] %d / %d probe cells pass CA1 depth filter\n', ...
                sum(ca1Mask), length(probeIdx));

        S_probe = subsetSpikes(S, ca1Idx);

        % Process all CA1 cells on this probe against this SWR file
        res = processCells(S_probe, probeSWR, session_duration, sessionID, ...
                           odor_to_outcome, SWR_PARTICIPATION_THRESH);

        % Append to accumulators
        all_session           = [all_session;           res.session];
        all_label             = [all_label;             res.label];
        all_location          = [all_location;          res.location];
        all_baseline_fr       = [all_baseline_fr;       res.baseline_fr];
        all_mean_fr_allSWR    = [all_mean_fr_allSWR;    res.mean_fr_allSWR];
        all_mean_fr_partSWR   = [all_mean_fr_partSWR;   res.mean_fr_partSWR];
        all_pct_participation = [all_pct_participation; res.pct_participation];
        all_participates      = [all_participates;      res.participates];
        all_swr_count         = [all_swr_count;         res.swr_count];
        all_odor              = [all_odor;              res.odor];
        all_odor_dir          = [all_odor_dir;          res.odor_dir];
        all_odor_outcome      = [all_odor_outcome;      res.odor_outcome];
        all_reward_dir        = [all_reward_dir;        res.reward_dir];
        all_airpuff_dir       = [all_airpuff_dir;       res.airpuff_dir];
        all_cell_type         = [all_cell_type;         res.cell_type];

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
    all_odor,              ...
    all_odor_dir,          ...
    all_odor_outcome,      ...
    all_reward_dir,        ...
    all_airpuff_dir,       ...
    all_baseline_fr,       ...
    all_mean_fr_allSWR,    ...
    all_mean_fr_partSWR,   ...
    all_pct_participation, ...
    all_participates,      ...
    all_swr_count,         ...
    'VariableNames', { ...
        'Session', 'Label', 'Location', 'CellType', ...
        'Odor', 'OdorDirection', 'OdorOutcome', ...
        'RewardDirection', 'AirpuffDirection', ...
        'BaselineFiringRate', 'MeanFiringRate_AllSWR', 'MeanFiringRate_ParticipatedSWR', ...
        'PctSWRParticipation', 'ParticipatesInSWR', 'SWRCount'});

filename     = char(strcat(ExpKeys.subject, '_', ExpKeys.date, '_SWR_participation_allcells.mat'));
filename_csv = char(strcat(ExpKeys.subject, '_', ExpKeys.date, '_SWR_participation_allcells.csv'));

save(filename, 'celltable', 'ExpKeys');
writetable(celltable, filename_csv);

fprintf('\nResults saved to:\n  %s\n  %s\n', filename, filename_csv);


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function csc = loadLFP(cfg, cfg_md, ExpKeys)
% loadLFP  Load the LFP for the requested probe, convert to mvdmlab tsd
%          format, apply true-depth correction, and sort channels dorsal->ventral.
%
%   cfg.probe  : 1 = imec0, 2 = imec1
%   cfg_md     : struct with field maxdepth (µm)
%   ExpKeys    : experiment keys struct

    csc = tsd;

    if cfg.probe == 1
        load('imec0_clean_lfp.mat');                          % loads imec0
        csc.tvec           = imec0.lfp_tvec';
        csc.data           = double(imec0.lfp_traces');
        csc.label          = cellstr(string(imec0.channel_ids).strip().extractAfter('ap#')');
        csc.cfg.depths     = imec0.depths;
        csc.cfg.shank_ids  = imec0.shank_ids;
        Fs = struct;
        Fs.SamplingFrequency   = double(imec0.lfp_fs);
        csc.cfg.hdr{1,1}   = Fs;
        % Convert to true depth
        csc.cfg.depths     = (cfg_md.maxdepth - csc.cfg.depths) ...
                              * cos(deg2rad(ExpKeys.probe1_insertion_roll));
    else
        load('imec1_clean_lfp.mat');                          % loads imec1
        csc.tvec           = imec1.lfp_tvec';
        csc.data           = double(imec1.lfp_traces');
        csc.label          = cellstr(string(imec1.channel_ids).strip().extractAfter('ap#')');
        csc.cfg.depths     = imec1.depths;
        csc.cfg.shank_ids  = imec1.shank_ids;
        Fs = struct;
        Fs.SamplingFrequency   = double(imec1.lfp_fs);
        csc.cfg.hdr{1,1}   = Fs;
        % Convert to true depth (probe 2 uses probe2_Depth)
        maxdepth2          = 1000 * ExpKeys.probe2_Depth;
        csc.cfg.depths     = (maxdepth2 - csc.cfg.depths) ...
                              * cos(deg2rad(ExpKeys.probe2_insertion_roll));
    end

    % Sort channels from most dorsal (smallest depth) to most ventral
    [~, sort_idx]      = sort(csc.cfg.depths);
    csc.data           = csc.data(sort_idx, :);
    csc.label          = csc.label(sort_idx);
    csc.cfg.depths     = csc.cfg.depths(sort_idx);
    csc.cfg.shank_ids  = csc.cfg.shank_ids(sort_idx);
    clear sort_idx;
end

function lfp = selectSWRLfp(csc, ExpKeys, cfg)
% selectSWRLfp  Return the SWR LFP (and optionally control LFP) for the
%               chosen probe region.
%
%   cfg.region     : 'dorsal' or 'ventral'  (case-insensitive)
%   cfg.useControl : true/false – also extract the control channel

    if strcmpi(cfg.region, 'dorsal')
        swrCh  = ExpKeys.imec0_dorsal_SWR_channel;
        ctrlCh = ExpKeys.imec0_dorsal_control_channel;
    else
        swrCh  = ExpKeys.imec0_ventral_SWR_channel;
        ctrlCh = ExpKeys.imec0_ventral_control_channel;
    end

    swrIdx         = find(cellfun(@(x) strcmp(x, swrCh),  csc.label), 1);
    lfp.SWR        = csc;
    lfp.SWR.data   = lfp.SWR.data(swrIdx, :);
    lfp.SWR.cfg.depths = lfp.SWR.cfg.depths(swrIdx);
    lfp.SWR.label  = lfp.SWR.label(swrIdx);

    if cfg.useControl
        ctrlIdx            = find(cellfun(@(x) strcmp(x, ctrlCh), csc.label), 1);
        lfp.Control        = csc;
        lfp.Control.data   = lfp.Control.data(ctrlIdx, :);
        lfp.Control.cfg.depths = lfp.Control.cfg.depths(ctrlIdx);
        lfp.Control.label  = lfp.Control.label(ctrlIdx);
    else
        lfp.Control = [];
    end
end

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
        % Fields that are not per-cell (e.g. scalars) are copied as-is
    end
end

function res = processCells(S, probeSWR, session_duration, sessionID, ...
                             odor_to_outcome, participation_thresh)
% processCells  Compute SWR participation metrics for every cell in S.
%
% All cells are processed regardless of selectivity. Cells with no
% selectivity information are labelled 'Unclassified'.

    nSWR = length(probeSWR.iv.tstart);

    res.session           = {};
    res.label             = {};
    res.location          = {};
    res.baseline_fr       = [];
    res.mean_fr_allSWR    = [];
    res.mean_fr_partSWR   = [];
    res.pct_participation = [];
    res.participates      = [];
    res.swr_count         = [];
    res.odor              = {};
    res.odor_dir          = {};
    res.odor_outcome      = {};
    res.reward_dir        = {};
    res.airpuff_dir       = {};
    res.cell_type         = {};

    for iCell = 1:length(S.label)

        % ---- Read selectivity fields ----------------------------------------
        odor_val        = S.usr.odor{iCell};
        odor_dir_val    = S.usr.odor_dir{iCell};
        reward_dir_val  = S.usr.reward_dir{iCell};
        airpuff_dir_val = S.usr.airpuff_dir{iCell};

        is_odor_selective     = ~strcmpi(odor_val,        'None') && ~isempty(odor_val);
        is_reward_responsive  = ~strcmpi(reward_dir_val,  'None') && ~isempty(reward_dir_val);
        is_airpuff_responsive = ~strcmpi(airpuff_dir_val, 'None') && ~isempty(airpuff_dir_val);
        is_outcome_responsive = is_reward_responsive || is_airpuff_responsive;

        % ---- Cell-type label -----------------------------------------------
        if is_odor_selective && is_outcome_responsive
            cell_type = 'Odor_and_Outcome';
        elseif is_odor_selective
            cell_type = 'Odor_Selective';
        elseif is_outcome_responsive
            cell_type = 'Outcome_Responsive';
        else
            cell_type = 'Unclassified';
        end

        % ---- Odor-associated outcome ----------------------------------------
        if is_odor_selective && isKey(odor_to_outcome, odor_val)
            odor_outcome_val = odor_to_outcome(odor_val);
        else
            odor_outcome_val = 'None';
        end

        % ---- SWR participation metrics --------------------------------------
        currCell_t               = S.t{iCell};
        swr_fr_matrix            = zeros(1, nSWR);
        swr_participation_binary = zeros(1, nSWR);

        for iSWR = 1:nSWR
            t_start = probeSWR.iv.tstart(iSWR);
            t_end   = probeSWR.iv.tend(iSWR);

            spike_count         = sum(currCell_t > t_start & currCell_t < t_end);
            swr_duration        = t_end - t_start;
            swr_fr_matrix(iSWR) = spike_count / swr_duration;

            if spike_count > 0
                swr_participation_binary(iSWR) = 1;
            end
        end

        mean_fr_allSWR = mean(swr_fr_matrix);

        part_idx = swr_participation_binary == 1;
        if any(part_idx)
            mean_fr_partSWR = mean(swr_fr_matrix(part_idx));
        else
            mean_fr_partSWR = 0;
        end

        pct_participation = mean(swr_participation_binary) * 100;
        participates      = pct_participation > participation_thresh;
        baseline_fr       = length(currCell_t) / session_duration;

        % ---- Append to result struct ----------------------------------------
        res.session           = [res.session;          {sessionID}];
        res.label             = [res.label;             {S.label{iCell}}];
        res.location          = [res.location;          {S.usr.HC_location{iCell}}];
        res.baseline_fr       = [res.baseline_fr;       baseline_fr];
        res.mean_fr_allSWR    = [res.mean_fr_allSWR;    mean_fr_allSWR];
        res.mean_fr_partSWR   = [res.mean_fr_partSWR;   mean_fr_partSWR];
        res.pct_participation = [res.pct_participation; pct_participation];
        res.participates      = [res.participates;      participates];
        res.swr_count         = [res.swr_count;         nSWR];
        res.odor              = [res.odor;              {odor_val}];
        res.odor_dir          = [res.odor_dir;          {odor_dir_val}];
        res.odor_outcome      = [res.odor_outcome;      {odor_outcome_val}];
        res.reward_dir        = [res.reward_dir;        {reward_dir_val}];
        res.airpuff_dir       = [res.airpuff_dir;       {airpuff_dir_val}];
        res.cell_type         = [res.cell_type;         {cell_type}];

    end % cell loop
end % processCells