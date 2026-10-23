%% Sharp Wave Ripple Participation Script - Dorsal vs Ventral Hippocampus
% Compares SWR participation by hippocampal location for only odor-selective
% and/or outcome-responsive cells.
% Selectivity/responsiveness read from new S_out.usr field structure:
%   usr.odor        - odor cell is selective for, or 'None'
%   usr.odor_dir    - 'excitatory', 'inhibitory', or 'None'
%   usr.reward_dir  - 'increase', 'decrease', or 'None'
%   usr.airpuff_dir - 'increase', 'decrease', or 'None'
%
% ADAPTED TO BE ABLE TO ANALYZE PARTICIPATION OF NEURONS IN SWRS THAT ARE
% ACTUALLY IN THEIR REGION (DORSAL vs. VENTRAL)
% EG May 2026 - adapted from odor_selective_swr_participation_dvHC3.m

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

SWR_MIN_COUNT = 20;        % Minimum number of SWRs required to include a session
SWR_PARTICIPATION_THRESH = 10; % % of SWRs a cell must fire in to count as "participating"

spikeFileBase = 'selective_permissive';

cfg_iv     = [];
cfg_iv.d   = [0, 0];  % Window (s) to expand each SWR interval by on each side

%% -------------------------------------------------------------------------
%% Initialize accumulators
%% -------------------------------------------------------------------------

all_session        = {};
all_label          = {};
all_location       = {};
all_baseline_fr    = [];
all_mean_fr_allSWR = [];
all_mean_fr_partSWR= [];
all_pct_participation = [];
all_participates   = [];
all_swr_count      = [];
all_odor           = {};
all_odor_dir       = {};
all_odor_outcome   = {};   % outcome associated with the preferred odor
all_reward_dir     = {};
all_airpuff_dir    = {};
all_cell_type      = {};
all_ExpKeys        = {};   % one ExpKeys struct per folder

%% -------------------------------------------------------------------------
%% Main loop - folders
%% -------------------------------------------------------------------------

for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});

    LoadExpKeys;
    all_ExpKeys{end+1} = ExpKeys;   % capture before any continue statements can skip it
    load all_events.mat;  % loads evt struct

    % ---- Depth calibration --------------------------------------------------
    record_depth   = ExpKeys.probe1_Depth;
    insertion_roll = ExpKeys.probe1_insertion_roll;
    cfg_md.maxdepth = 1000 * record_depth;

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

    % ---- Find all SWR files in this folder ----------------------------------
    AllSWR = FindFiles('*SWR.mat');

    % loops through the SWR files for a session
    for currFile = 1:length(AllSWR)

        swrFile = AllSWR{currFile};

        % Determine probe ID
        if contains(swrFile, 'imec0')
            probeID   = 'imec0';
            % spikefile = 'selective_imec0_all.mat';
            spikefile = strcat(spikeFileBase, '_', probeID, '_all.mat');
            depth_offset = 0;
            DH_threshold = 3000;
        elseif contains(swrFile, 'imec1')
            probeID   = 'imec1';
            % spikefile = 'selective_imec1_all.mat';
            spikefile = strcat(spikeFileBase, '_', probeID, '_all.mat');
            depth_offset = 0;   % additional offset for imec1
            DH_threshold = 3000;
        else
            continue;  % skip files that don't match a known probe
        end

        % Check spike data exists
        if exist(spikefile, 'file') ~= 2
            fprintf('  [skip] %s not found\n', spikefile);
            continue;
        end

        % Load SWR intervals
        load(swrFile);
        if contains(swrFile, 'dorsal')
            probeSWR = dorsal_swr;
            swr_location = 'Dorsal';
        elseif contains(swrFile, 'ventral')
            probeSWR = ventral_swr;
            swr_location = 'Ventral';
        else
            fprintf('  [skip] Cannot determine dorsal/ventral from filename: %s\n', swrFile);
            continue;
        end
        probeSWR.probeID = probeID;
        probeSWR.iv      = ResizeIV(cfg_iv, probeSWR.iv);

        % Skip sessions with too few SWRs
        if length(probeSWR.iv.tstart) < SWR_MIN_COUNT
            fprintf('  [skip] Too few SWRs (%d) in %s\n', length(probeSWR.iv.tstart), swrFile);
            continue;
        end

        % Load spike data
        load(spikefile);   % loads S_all
        S = S_out;

        % Compute true depth from probe tip and classify dorsal vs ventral
        S.usr.true_depth = (cfg_md.maxdepth - S.usr.recording_depths - depth_offset) ...
                            .* cos(deg2rad(insertion_roll));
        for iCell = 1:length(S.label)
            if S.usr.true_depth(iCell) < DH_threshold
                S.usr.HC_location{iCell} = 'Dorsal';
            else
                S.usr.HC_location{iCell} = 'Ventral';
            end
        end

        % Process all cells for this probe/session
        sessionID = strcat(ExpKeys.subject, '-', ExpKeys.date);
        [res] = processCells(S, probeSWR, session_duration, sessionID, ...
                             odor_to_outcome, SWR_PARTICIPATION_THRESH, ...
                             swr_location);

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

if isscalar(folder_list)
    filename     = char(strcat(ExpKeys.subject, '_', ExpKeys.date, '_SWR_participation.mat'));
    filename_csv = char(strcat(ExpKeys.subject, '_', ExpKeys.date, '_SWR_participation.csv'));
elseif length(folder_list) > 1
    folderPath = 'D:\vCA1';
    baseName = 'allSessions_SWR_participation';
    filename     = fullfile(folderPath, strcat(baseName, '.mat'));
    filename_csv = fullfile(folderPath, strcat(baseName, '.csv'));
end

save(filename, 'celltable', 'all_ExpKeys', 'cfg_iv');
writetable(celltable, filename_csv);

fprintf('\nResults saved to:\n  %s\n  %s\n', filename, filename_csv);


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function res = processCells(S, probeSWR, session_duration, sessionID, ...
                             odor_to_outcome, participation_thresh, swr_location)
% processCells  Compute SWR participation metrics for every qualifying cell
% in spike struct S against the SWR intervals in probeSWR.
%
% Inputs:
%   S                    - spike struct (S_all), with usr fields:
%                            .odor, .odor_dir, .reward_dir, .airpuff_dir,
%                            .HC_location, .recording_depths, .true_depth
%   probeSWR             - SWR struct with .iv.tstart / .iv.tend
%   session_duration     - scalar, length of session in seconds
%   sessionID            - string, e.g. 'M649-2026_03_30'
%   odor_to_outcome      - containers.Map: odor name -> outcome string
%   participation_thresh - scalar, % threshold for binary participation
%   swr_location         - string, 'Dorsal' or 'Ventral'; only cells from
%                            the matching HC region are processed
%
% Output:
%   res - struct of column vectors/cell arrays, one entry per qualifying cell

    nSWR = length(probeSWR.iv.tstart);

    % Pre-allocate result struct as empty cell/numeric arrays
    res.session        = {};
    res.label          = {};
    res.location       = {};
    res.baseline_fr    = [];
    res.mean_fr_allSWR = [];
    res.mean_fr_partSWR= [];
    res.pct_participation = [];
    res.participates   = [];
    res.swr_count      = [];
    res.odor           = {};
    res.odor_dir       = {};
    res.odor_outcome   = {};
    res.reward_dir     = {};
    res.airpuff_dir    = {};
    res.cell_type      = {};

    for iCell = 1:length(S.label)

        % ---- Skip cells not in the matching HC region ----------------------
        if ~strcmpi(S.usr.HC_location{iCell}, swr_location)
            continue;
        end

        % ---- Read selectivity / responsiveness fields -----------------------
        odor_val        = S.usr.odor{iCell};
        odor_dir_val    = S.usr.odor_dir{iCell};
        reward_dir_val  = S.usr.reward_dir{iCell};
        airpuff_dir_val = S.usr.airpuff_dir{iCell};

        is_odor_selective   = ~strcmpi(odor_val,        'None') && ~isempty(odor_val);
        is_reward_responsive  = ~strcmpi(reward_dir_val,  'None') && ~isempty(reward_dir_val);
        is_airpuff_responsive = ~strcmpi(airpuff_dir_val, 'None') && ~isempty(airpuff_dir_val);
        is_outcome_responsive = is_reward_responsive || is_airpuff_responsive;

        % Skip cells with no relevant selectivity
        if ~is_odor_selective && ~is_outcome_responsive
            continue;
        end

        % ---- Determine cell-type label -------------------------------------
        if is_odor_selective && is_outcome_responsive
            cell_type = 'Odor_and_Outcome';
        elseif is_odor_selective
            cell_type = 'Odor_Selective';
        else
            cell_type = 'Outcome_Responsive';
        end

        % ---- Odor-associated outcome ---------------------------------------
        if is_odor_selective && isKey(odor_to_outcome, odor_val)
            odor_outcome_val = odor_to_outcome(odor_val);
        else
            odor_outcome_val = 'None';
        end

        % ---- SWR participation metrics -------------------------------------
        currCell_t          = S.t{iCell};
        swr_fr_matrix       = zeros(1, nSWR);
        swr_participation_binary = zeros(1, nSWR);

        for iSWR = 1:nSWR
            t_start = probeSWR.iv.tstart(iSWR);
            t_end   = probeSWR.iv.tend(iSWR);

            spike_count  = sum(currCell_t > t_start & currCell_t < t_end);
            swr_duration = t_end - t_start;
            swr_fr_matrix(iSWR) = spike_count / swr_duration;

            if spike_count > 0
                swr_participation_binary(iSWR) = 1;
            end
        end

        % Mean FR across ALL SWRs (including non-participation SWRs)
        mean_fr_allSWR    = mean(swr_fr_matrix);

        % Mean FR restricted to SWRs the cell actually fired in
        part_idx = swr_participation_binary == 1;
        if any(part_idx)
            mean_fr_partSWR = mean(swr_fr_matrix(part_idx));
        else
            mean_fr_partSWR = 0;
        end

        % Percent of SWRs with at least one spike
        pct_participation = mean(swr_participation_binary) * 100;

        % Binary: does the cell clear the participation threshold?
        participates = pct_participation > participation_thresh;

        % Baseline FR over full session
        baseline_fr = length(currCell_t) / session_duration;

        % ---- Append to result struct ---------------------------------------
        res.session           = [res.session;        {sessionID}];
        res.label             = [res.label;           {S.label{iCell}}];
        res.location          = [res.location;        {S.usr.HC_location{iCell}}];
        res.baseline_fr       = [res.baseline_fr;     baseline_fr];
        res.mean_fr_allSWR    = [res.mean_fr_allSWR;  mean_fr_allSWR];
        res.mean_fr_partSWR   = [res.mean_fr_partSWR; mean_fr_partSWR];
        res.pct_participation = [res.pct_participation; pct_participation];
        res.participates      = [res.participates;    participates];
        res.swr_count         = [res.swr_count;       nSWR];
        res.odor              = [res.odor;            {odor_val}];
        res.odor_dir          = [res.odor_dir;        {odor_dir_val}];
        res.odor_outcome      = [res.odor_outcome;    {odor_outcome_val}];
        res.reward_dir        = [res.reward_dir;      {reward_dir_val}];
        res.airpuff_dir       = [res.airpuff_dir;     {airpuff_dir_val}];
        res.cell_type         = [res.cell_type;       {cell_type}];

    end % cell loop
end % processCells