%% SWR Analysis Script — multi-session wrapper 
% Analyzes sharp-wave ripple (SWR) occurrence between experimental trials.
% Built off of swr_experience_dependence_v5.m
% Normalization:
%   For each session and each probe, every ITI's raw SWR count is divided by
%   the session-wide total SWR count for that probe before dividing by ITI
%   duration, then multiplied by 100.  The resulting "normalized rate" sums
%   to 100 across all ITIs of a session, so bars show each trial-type's
%   percentage share of the session's total SWRs rather than absolute
%   events/s. For between session analysis, it averages the average. So
%   basically, because there aren't necessarily an equal number of ITIs for
%   each trial type, the combined average is the average of the average
%   normalized rate for each trial type in a single session.
% 
% Session mask: Doesn't analyze regions with fewer than 60 SWRs in that
% region during a single session
%
% Experiment structure:
%   - Sequential odor pairs followed by outcome (reward, airpuff, or neutral)
%   - Possible catch trials (outcome without odors, or odors without outcome)
%
% Trial types classified:
%   1. reward              - reward odors + reward outcome
%   2. airpuff             - airpuff odors + airpuff outcome
%   3. neutral             - neutral odors + neutral outcome
%   4. reward-catch        - reward outcome with NO preceding odors
%   5. airpuff-catch       - airpuff outcome with NO preceding odors
%   6. odor-reward-catch   - reward odors with NO subsequent reward
%   7. odor-airpuff-catch  - airpuff odors with NO subsequent airpuff
%
% Outputs (saved automatically):
%   <session_dir>/SWR_ITI_results_minITI_trialSpecific.mat       — per-session ITI table + summary
%   <session_dir>/SWR_ITI_plot_minITI_trialSpecific.pdf          — per-session bar figure (raw)
%   <session_dir>/SWR_ITI_plot_norm_minITI_trialSpecific.pdf     — per-session bar figure (normalized)
%   <output_root>/SWR_combined_results_minITI_trialSpecific.mat  — pooled data across sessions
%   <output_root>/SWR_combined_plot_minITI_trialSpecific.pdf     — pooled bar figure (raw)
%   <output_root>/SWR_combined_plot_norm_minITI_trialSpecific.pdf— pooled bar figure (normalized)

%% ========================================================================
%  CONFIGURATION — edit these before running
%  ========================================================================
% Just good unit sessions
% session_dirs = {
%     'D:\vCA1\M560\preprocessed\M560_2025_12_19'
%     'D:\vCA1\M593\preprocessed\M593_2025_12_26'
%     'D:\vCA1\M593\preprocessed\M593_2025_12_28'
%     'D:\vCA1\M653\preprocessed\M653_2026_02_04'
%     'D:\vCA1\M656\preprocessed\M656_2026_03_25'
%     'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
%     'D:\vCA1\M649\preprocessed\M649_2026_04_01'
%     % add more sessions here ...
% };

% All SWR sessions
session_dirs = {
    'D:\vCA1\M560\preprocessed\M560_2025_12_19'
    'D:\vCA1\M593\preprocessed\M593_2025_12_26'
    'D:\vCA1\M593\preprocessed\M593_2025_12_28'
    'D:\vCA1\M653\preprocessed\M653_2026_02_04'
    'D:\vCA1\M656\preprocessed\M656_2026_03_25'
    'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
    'D:\vCA1\M649\preprocessed\M649_2026_04_01'
    'D:\vCA1\M653\preprocessed\M653_2026_02_05'
    'D:\vCA1\M653\preprocessed\M653_2026_02_07'
    'D:\vCA1\M656\preprocessed\M656_2026_03_23'
    'D:\vCA1\M656\preprocessed\M656_2026_03_24'
    'D:\vCA1\M649\preprocessed\M649_2026_04_02'
    'D:\vCA1\M649\preprocessed\M649_2026_04_03'
    % add more sessions here ...
};

% session_dirs = {
%     % add more sessions here ...
%     'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
%     'D:\vCA1\M649\preprocessed\M649_2026_04_01'
%     'D:\vCA1\M649\preprocessed\M649_2026_04_02'
%     'D:\vCA1\M649\preprocessed\M649_2026_04_03'
% };

% output_root = 'D:\vCA1\M593\preprocessed\combined_results_normalized';
% output_root = 'D:\vCA1\M653\preprocessed\M653_2026_02_05\SWR_results';
output_root = 'D:\vCA1\graveyard';


plotSession = false; % toggle to true or false depending on whether you want all the session plots to be output

MIN_ITI_DURATION     = 5.0;   % (s) ITIs shorter than this are excluded
MIN_SWR_COUNT        = 60;    % minimum total SWRs in a session for that region to be included. Set to 0 if you don't want

ODOR_TO_OUTCOME_MIN  = 2.7;   % (s) odor 1 onset to outcome onset
ODOR_TO_OUTCOME_MAX  = 3.05;  % (s) odor 1 onset to outcome onset
CATCH_TRIAL_WINDOW   = 2.5;   % (s) how long after odor 1 onset an odor-catch trial ends (second odor offset)
NEUTRAL_TRIAL_WINDOW = 2.5;   % (s) how long after odor 1 onset a neutral trial ends (second odor offset)
REWARD_TRIAL_BUFFER  = 10;    % (s) buffer after reward admin before ITI starts
AIRPUFF_TRIAL_BUFFER = 5;     % (s) buffer after airpuff admin before ITI starts

addpath(genpath('C:\Users\mvdmlab\Documents\egibbens\odor-pixels\analysis\matlab'));
addpath(genpath('C:\Users\mvdmlab\Documents\egibbens\EG_analysis\odor-pixels\analysis\matlab'));

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

%% ========================================================================
%  SECTION 1: Iterate over sessions
%  ========================================================================
all_sessions   = struct();
combined_ITI   = struct([]);
combined_count = 0;

canonical_types = {'reward','airpuff','neutral', ...
                   'reward-catch','airpuff-catch', ...
                   'odor-reward-catch','odor-airpuff-catch'};

for sIdx = 1:numel(session_dirs)

    sess_dir = session_dirs{sIdx};
    [~, sess_name] = fileparts(sess_dir);
    fprintf('\n\n================================================================\n');
    fprintf('  Processing session %d / %d : %s\n', sIdx, numel(session_dirs), sess_name);
    fprintf('================================================================\n');

    %% --------------------------------------------------------------------
    %  SECTION 1a: Load session data
    %  --------------------------------------------------------------------
    cd(sess_dir);
    LoadExpKeys;

    load('all_events.mat');
    fprintf('Loaded all_events.mat\n');

    % --- Lost-signal epochs (Nx2 matrix of [start, end] pairs in seconds).
    %     Defaults to empty so sessions without this field are unaffected.
    if isfield(ExpKeys, 'imec0_lost_signal_times')
        lost_signal_times = ExpKeys.imec0_lost_signal_times; % Nx2
        fprintf('Lost-signal epochs loaded: %d period(s).\n', size(lost_signal_times, 1));
    else
        lost_signal_times = zeros(0, 2); % empty — no exclusions
    end

    has_dorsal  = false;
    has_ventral = false;

    if isfile('dorsal_imec0_SWR.mat')
        load('dorsal_imec0_SWR.mat');
        has_dorsal = true;
        fprintf('Loaded dorsal_imec0_SWR.mat\n');
    else
        warning('Session %s: dorsal_imec0_SWR.mat not found — dorsal SWRs skipped.', sess_name);
    end

    if isfile('ventral_imec0_SWR.mat')
        load('ventral_imec0_SWR.mat');
        has_ventral = true;
        fprintf('Loaded ventral_imec0_SWR.mat\n');
    else
        warning('Session %s: ventral_imec0_SWR.mat not found — ventral SWRs skipped.', sess_name);
    end

    if ~has_dorsal && ~has_ventral
        warning('Session %s: no SWR files found — skipping session.', sess_name);
        continue
    end

    %% --------------------------------------------------------------------
    %  SECTION 1b: Extract SWR event times
    %  --------------------------------------------------------------------
    if has_dorsal
        [dTstart, dTend] = extractSWRtimes(dorsal_swr, 'dorsal_swr');
        dorsal_SWR_times = (dTstart + dTend) / 2;
        fprintf('Dorsal  SWR count: %d\n', numel(dorsal_SWR_times));

        % reject dorsal if total SWR count is below threshold
        if numel(dorsal_SWR_times) < MIN_SWR_COUNT
            warning('Session %s: dorsal SWR count (%d) below MIN_SWR_COUNT (%d) — dorsal SWRs rejected.', ...
                sess_name, numel(dorsal_SWR_times), MIN_SWR_COUNT);
            has_dorsal       = false;
            dorsal_SWR_times = [];
        end
    else
        dorsal_SWR_times = [];
    end

    if has_ventral
        [vTstart, vTend]  = extractSWRtimes(ventral_swr, 'ventral_swr');
        ventral_SWR_times = (vTstart + vTend) / 2;
        fprintf('Ventral SWR count: %d\n', numel(ventral_SWR_times));

        % reject ventral if total SWR count is below threshold
        if numel(ventral_SWR_times) < MIN_SWR_COUNT
            warning('Session %s: ventral SWR count (%d) below MIN_SWR_COUNT (%d) — ventral SWRs rejected.', ...
                sess_name, numel(ventral_SWR_times), MIN_SWR_COUNT);
            has_ventral       = false;
            ventral_SWR_times = [];
        end
    else
        ventral_SWR_times = [];
    end

    % skip session if both regions were rejected
    if ~has_dorsal && ~has_ventral
        warning('Session %s: both regions rejected due to low SWR count — skipping session.', sess_name);
        continue
    end

    %% --------------------------------------------------------------------
    %  SECTION 1c: Parse odor / outcome events
    %  --------------------------------------------------------------------
    reward_odors  = ExpKeys.reward_odors;
    airpuff_odors = ExpKeys.airpuff_odors;
    neutral_odors = ExpKeys.neutral_odors;

    reward_odor1_times  = getEventTimes(evt, reward_odors{1});
    airpuff_odor1_times = getEventTimes(evt, airpuff_odors{1});
    neutral_odor1_times = getEventTimes(evt, neutral_odors{1});

    reward_times  = getEventTimes(evt, 'Reward ON');
    airpuff_times = getEventTimes(evt, 'Airpuff ON');

    %% --------------------------------------------------------------------
    %  SECTION 1d: Classify trials
    %  --------------------------------------------------------------------
    trials = classifyTrials( ...
        reward_odor1_times,  airpuff_odor1_times, neutral_odor1_times, ...
        reward_times,        airpuff_times, ...
        ODOR_TO_OUTCOME_MIN, ODOR_TO_OUTCOME_MAX, ...
        CATCH_TRIAL_WINDOW,  NEUTRAL_TRIAL_WINDOW, REWARD_TRIAL_BUFFER, AIRPUFF_TRIAL_BUFFER);

    [~, sortIdx] = sort([trials.startTime]);
    trials = trials(sortIdx);

    trialTypes = {trials.type};
    uTypes     = unique(trialTypes);

    fprintf('\n=== Trial classification summary ===\n');
    for i = 1:numel(uTypes)
        fprintf('  %-25s : %d trials\n', uTypes{i}, sum(strcmp(trialTypes, uTypes{i})));
    end

    %% --------------------------------------------------------------------
    %  SECTION 1e: Compute raw SWR rate in each ITI
    %  --------------------------------------------------------------------
    nTrials     = numel(trials);
    ITI_results = struct();
    itiCount    = 0;
    nExcluded   = 0;
    nExcludedSignal = 0; % ITIs excluded due to lost-signal overlap

    for i = 1:(nTrials - 1)

        t_start  = trials(i).endTime;
        t_end    = trials(i+1).startTime;
        duration = t_end - t_start;

        if duration < MIN_ITI_DURATION
            nExcluded = nExcluded + 1;
            continue
        end
        
        % --- Skip ITIs that overlap with any lost-signal epoch.
        %     Overlap condition: ITI starts before epoch ends AND
        %                        ITI ends after epoch starts.
        if ~isempty(lost_signal_times)
            overlaps = t_start < lost_signal_times(:,2) & ...
                       t_end   > lost_signal_times(:,1);
            if any(overlaps)
                nExcludedSignal = nExcludedSignal + 1;
                continue
            end
        end

        itiCount = itiCount + 1;

        if has_dorsal
            n_dorsal    = sum(dorsal_SWR_times  >= t_start & dorsal_SWR_times  <= t_end);
            rate_dorsal = n_dorsal / duration;
        else
            n_dorsal    = NaN;
            rate_dorsal = NaN;
        end

        if has_ventral
            n_ventral    = sum(ventral_SWR_times >= t_start & ventral_SWR_times <= t_end);
            rate_ventral = n_ventral / duration;
        else
            n_ventral    = NaN;
            rate_ventral = NaN;
        end

        ITI_results(itiCount).session           = sess_name;
        ITI_results(itiCount).trialIndex_before = i;
        ITI_results(itiCount).trialIndex_after  = i + 1;
        ITI_results(itiCount).type_before       = trials(i).type;
        ITI_results(itiCount).type_after        = trials(i+1).type;
        ITI_results(itiCount).t_start           = t_start;
        ITI_results(itiCount).t_end             = t_end;
        ITI_results(itiCount).duration_s        = duration;
        ITI_results(itiCount).n_dorsal_SWR      = n_dorsal;
        ITI_results(itiCount).n_ventral_SWR     = n_ventral;
        ITI_results(itiCount).rate_dorsal_SWR   = rate_dorsal;
        ITI_results(itiCount).rate_ventral_SWR  = rate_ventral;
        % Normalized fields are filled in Section 1e-norm below
        ITI_results(itiCount).rate_dorsal_SWR_norm  = NaN;
        ITI_results(itiCount).rate_ventral_SWR_norm = NaN;
    end

    fprintf('\n=== ITI analysis complete ===\n');
    fprintf('  ITIs included : %d\n', itiCount);
    fprintf('  ITIs excluded (< %.1f s) : %d\n', MIN_ITI_DURATION, nExcluded);
    fprintf('  ITIs excluded (lost signal)       : %d\n', nExcludedSignal);

    %% --------------------------------------------------------------------
    %  SECTION 1e-norm: Normalize rates by session-wide SWR total
    %
    %  For each region, sum all raw counts across every included ITI in this
    %  session.  Then divide each ITI's count by that session total and by
    %  the ITI duration, giving a normalized rate whose integral over all
    %  ITIs equals 100 % of the session's SWRs.
    %
    %  If the session total is 0 (no SWRs detected at all for a region),
    %  normalized rates are left as NaN to avoid division-by-zero artefacts.
    %  --------------------------------------------------------------------
    if has_dorsal
        sess_total_dorsal = sum([ITI_results.n_dorsal_SWR], 'omitnan');
        fprintf('  Session dorsal  SWR total (in included ITIs): %d\n', sess_total_dorsal);
    end
    if has_ventral
        sess_total_ventral = sum([ITI_results.n_ventral_SWR], 'omitnan');
        fprintf('  Session ventral SWR total (in included ITIs): %d\n', sess_total_ventral);
    end

    for i = 1:itiCount
        dur = ITI_results(i).duration_s;

        if has_dorsal && sess_total_dorsal > 0
            ITI_results(i).rate_dorsal_SWR_norm = ...
                100 * ITI_results(i).n_dorsal_SWR / (sess_total_dorsal * dur);
        end

        if has_ventral && sess_total_ventral > 0
            ITI_results(i).rate_ventral_SWR_norm = ...
                100 * ITI_results(i).n_ventral_SWR / (sess_total_ventral * dur);
        end
    end

    %% --------------------------------------------------------------------
    %  SECTION 1f: Per-session summary statistics (printed + saved)
    %  --------------------------------------------------------------------
    sess_summary      = computeSummary(ITI_results, canonical_types, has_dorsal, has_ventral);
    sess_summary_norm = computeSummaryNorm(ITI_results, canonical_types, has_dorsal, has_ventral);

    fprintf('\n--- Per-session raw SWR rate (SWR/s) ---\n');
    printSummary(sess_summary, has_dorsal, has_ventral);

    fprintf('\n--- Per-session normalized SWR rate (percentage of session SWRs / s) ---\n');
    printSummaryNorm(sess_summary_norm, has_dorsal, has_ventral);

    save_path_mat = fullfile(sess_dir, 'SWR_ITI_results_minITI_trialSpecific.mat');
    save(save_path_mat, 'ITI_results', 'sess_summary', 'sess_summary_norm', ...
        'has_dorsal', 'has_ventral', 'MIN_ITI_DURATION', 'REWARD_TRIAL_BUFFER', ...
        'AIRPUFF_TRIAL_BUFFER', 'sess_name');
    fprintf('Saved: %s\n', save_path_mat);

    %% --------------------------------------------------------------------
    %  SECTION 1g: Per-session plots (raw + normalized)
    %  --------------------------------------------------------------------
    if plotSession
        fig_sess = plotSWRbars(ITI_results, canonical_types, has_dorsal, has_ventral, ...
            sprintf('Session: %s', strrep(sess_name,'_',' ')));
        save_path_fig = fullfile(sess_dir, 'SWR_ITI_plot_minITI_trialSpecific');
        exportgraphics(fig_sess, strcat(save_path_fig, '.pdf'), 'ContentType', 'vector');
        exportgraphics(fig_sess, strcat(save_path_fig, '.jpg'), 'Resolution', 300);
        fprintf('Saved: %s\n', save_path_fig);
        close(fig_sess);
    
        fig_sess_norm = plotSWRbarsNorm(ITI_results, canonical_types, has_dorsal, has_ventral, ...
            sprintf('Session: %s  (normalized)', strrep(sess_name,'_',' ')));
        save_path_fig_norm = fullfile(sess_dir, 'SWR_ITI_plot_norm_minITI_trialSpecific');
        exportgraphics(fig_sess_norm, strcat(save_path_fig_norm, '.pdf'), 'ContentType', 'vector');
        exportgraphics(fig_sess_norm, strcat(save_path_fig_norm, '.jpg'), 'Resolution', 300);
        fprintf('Saved: %s\n', save_path_fig_norm);
        close(fig_sess_norm);
    end

    %% --------------------------------------------------------------------
    %  SECTION 1h: Accumulate into combined pool
    %  --------------------------------------------------------------------
    if isempty(fieldnames(combined_ITI))
        combined_ITI   = ITI_results;
        combined_count = numel(ITI_results);
    else
        combined_ITI   = [combined_ITI, ITI_results];
        combined_count = numel(combined_ITI);
    end

    all_sessions(sIdx).name        = sess_name;
    all_sessions(sIdx).has_dorsal  = has_dorsal;
    all_sessions(sIdx).has_ventral = has_ventral;
    all_sessions(sIdx).nITI        = itiCount;
    all_sessions(sIdx).nExcluded   = nExcluded;
    all_sessions(sIdx).summary     = sess_summary;
    all_sessions(sIdx).summary_norm = sess_summary_norm;

    %% --------------------------------------------------------------------
    %  SECTION 1i: Clean up session-specific variables
    %  --------------------------------------------------------------------
    clearvars -except ...
        session_dirs output_root plotSession ...
        MIN_ITI_DURATION MIN_SWR_COUNT ...
        ODOR_TO_OUTCOME_MIN ODOR_TO_OUTCOME_MAX ...
        CATCH_TRIAL_WINDOW NEUTRAL_TRIAL_WINDOW REWARD_TRIAL_BUFFER ...
        AIRPUFF_TRIAL_BUFFER canonical_types ...
        all_sessions combined_ITI combined_count ...
        sIdx

end  % end session loop

%% ========================================================================
%  SECTION 2: Combined (all sessions) analysis
%  ========================================================================
if combined_count == 0
    error('No valid ITI data found across any session.');
end

any_dorsal  = any([all_sessions.has_dorsal]);
any_ventral = any([all_sessions.has_ventral]);

fprintf('\n\n================================================================\n');
fprintf('  COMBINED RESULTS  (%d sessions, %d ITIs)\n', numel(session_dirs), combined_count);
fprintf('================================================================\n');

combined_summary      = computeSummary(combined_ITI, canonical_types, any_dorsal, any_ventral);
combined_summary_norm = computeSummaryNormHierarchical(all_sessions, canonical_types, any_dorsal, any_ventral);

printSummary(combined_summary, any_dorsal, any_ventral);
printSummaryNorm(combined_summary_norm, any_dorsal, any_ventral);

save_path_combined = fullfile(output_root, 'SWR_combined_results_minITI_trialSpecific.mat');
save(save_path_combined, 'combined_ITI', 'combined_summary', 'combined_summary_norm', ...
    'all_sessions', 'canonical_types', 'MIN_ITI_DURATION', 'MIN_SWR_COUNT', ...
    'REWARD_TRIAL_BUFFER', 'AIRPUFF_TRIAL_BUFFER');
fprintf('\nSaved combined results: %s\n', save_path_combined);

% Raw combined plot — unchanged, ITI-level pooling is fine for raw rates
fig_combined = plotSWRbars(combined_ITI, canonical_types, any_dorsal, any_ventral, ...
    sprintf('All sessions combined  (n = %d ITIs)', combined_count));
save_path_combined_fig = fullfile(output_root, 'SWR_combined_plot_minITI_trialSpecific');
exportgraphics(fig_combined, strcat(save_path_combined_fig, '.pdf'), 'ContentType', 'vector');
exportgraphics(fig_combined, strcat(save_path_combined_fig, '.jpg'), 'Resolution', 300);
fprintf('Saved combined plot       : %s\n', save_path_combined_fig);
close(fig_combined);

% Normalized combined plot — now uses hierarchical summary
n_sessions = numel(all_sessions);
fig_combined_norm = plotSWRbarsNormHierarchical(combined_summary_norm, canonical_types, ...
    any_dorsal, any_ventral, ...
    sprintf('All sessions combined — normalized  (n = %d sessions)', n_sessions));
save_path_combined_fig_norm = fullfile(output_root, 'SWR_combined_plot_norm_minITI_trialSpecific');
exportgraphics(fig_combined_norm, strcat(save_path_combined_fig_norm, '.pdf'), 'ContentType', 'vector');
exportgraphics(fig_combined_norm, strcat(save_path_combined_fig_norm, '.jpg'), 'Resolution', 300);
fprintf('Saved combined norm plot  : %s\n', save_path_combined_fig_norm);
close(fig_combined_norm);

fprintf('\nAll done.\n');

%% ========================================================================
%  SECTION 3: Statistical analyses
%  
%
%  Unit of observation : one SESSION-LEVEL MEAN per (session × region ×
%                        trial-type) cell, using the hierarchical
%                        session-means already stored in combined_summary_norm.
%                        This is consistent with the bar-plot averaging
%                        strategy and avoids pseudo-replication from ITIs.
%
%  Only the three *core* trial types are included:
%       'neutral'  'reward'  'airpuff'
%
%  Animal ID is extracted as the first four characters of the session name
%  (e.g.  'M560_2025_12_19'  →  'M560').
%  ========================================================================

core_types = {'neutral', 'reward', 'airpuff'};   % Factor B levels
n_sessions = numel(all_sessions);

% -------------------------------------------------------------------------
%  SECTION 3a: Build the long-format data table
%  Columns: SWR_rate | Region | TrialType | Animal | Session
% -------------------------------------------------------------------------
rows_rate      = [];
rows_region    = {};
rows_trialtype = {};
rows_animal    = {};
rows_session   = {};

for s = 1:n_sessions
    sess_nm   = all_sessions(s).name;
    animal_id = sess_nm(1:4);          % first 4 chars  →  e.g. 'M560'
    sess_norm = all_sessions(s).summary_norm;

    for ti = 1:numel(core_types)
        ttype      = core_types{ti};
        type_match = strcmp({sess_norm.type}, ttype);
        if ~any(type_match), continue; end

        % --- Dorsal ---
        if all_sessions(s).has_dorsal
            dN   = sess_norm(type_match).dorsal_N;
            dMu  = sess_norm(type_match).dorsal_mean;
            if ~isnan(dN) && dN > 0
                rows_rate(end+1,1)      = dMu;        %#ok<SAGROW>
                rows_region{end+1,1}    = 'Dorsal';   %#ok<SAGROW>
                rows_trialtype{end+1,1} = ttype;      %#ok<SAGROW>
                rows_animal{end+1,1}    = animal_id;  %#ok<SAGROW>
                rows_session{end+1,1}   = sess_nm;    %#ok<SAGROW>
            end
        end

        % --- Ventral ---
        if all_sessions(s).has_ventral
            vN   = sess_norm(type_match).ventral_N;
            vMu  = sess_norm(type_match).ventral_mean;
            if ~isnan(vN) && vN > 0
                rows_rate(end+1,1)      = vMu;
                rows_region{end+1,1}    = 'Ventral';
                rows_trialtype{end+1,1} = ttype;
                rows_animal{end+1,1}    = animal_id;
                rows_session{end+1,1}   = sess_nm;
            end
        end
    end
end

% Convert to a MATLAB table with categorical predictors
stats_tbl = table( ...
    rows_rate,                           ...
    categorical(rows_region),            ...
    categorical(rows_trialtype),         ...
    categorical(rows_animal),            ...
    rows_session,                        ...
    'VariableNames', {'NormSWR_rate', 'Region', 'TrialType', 'Animal', 'Session'});

fprintf('\n\n================================================================\n');
fprintf('  SECTION 3 — Statistical analyses\n');
fprintf('  Long-format table: %d rows  (%d sessions × regions × trial types)\n', ...
    height(stats_tbl), height(stats_tbl));
fprintf('================================================================\n');
disp(stats_tbl);

% Save the table for external use
save(fullfile(output_root, 'SWR_stats_table_minITI_trialSpecific.mat'), 'stats_tbl');
writetable(stats_tbl, fullfile(output_root, 'SWR_stats_table_trialSpecific.csv'));
fprintf('Saved stats table: %s\n', fullfile(output_root, 'SWR_stats_table_minITI_trialSpecific.csv'));



%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

% -------------------------------------------------------------------------
function summary = computeSummary(ITI_results, types, has_dorsal, has_ventral)
% Returns a struct array (one entry per trial type) with N, mean, SEM,
% and total raw SWR count for dorsal and ventral probes.

summary = struct();
for i = 1:numel(types)
    summary(i).type = types{i};
    idx = strcmp({ITI_results.type_before}, types{i});

    if has_dorsal
        vals = [ITI_results(idx).rate_dorsal_SWR];
        n    = sum(~isnan(vals));
        summary(i).dorsal_N     = n;
        summary(i).dorsal_mean  = mean(vals, 'omitnan');
        summary(i).dorsal_sem   = std(vals, 'omitnan') / max(sqrt(n), 1);
        summary(i).dorsal_total = sum([ITI_results(idx).n_dorsal_SWR], 'omitnan');
    else
        summary(i).dorsal_N     = NaN;
        summary(i).dorsal_mean  = NaN;
        summary(i).dorsal_sem   = NaN;
        summary(i).dorsal_total = NaN;
    end

    if has_ventral
        vals = [ITI_results(idx).rate_ventral_SWR];
        n    = sum(~isnan(vals));
        summary(i).ventral_N     = n;
        summary(i).ventral_mean  = mean(vals, 'omitnan');
        summary(i).ventral_sem   = std(vals, 'omitnan') / max(sqrt(n), 1);
        summary(i).ventral_total = sum([ITI_results(idx).n_ventral_SWR], 'omitnan');
    else
        summary(i).ventral_N     = NaN;
        summary(i).ventral_mean  = NaN;
        summary(i).ventral_sem   = NaN;
        summary(i).ventral_total = NaN;
    end
end
end

% -------------------------------------------------------------------------
function summary = computeSummaryNorm(ITI_results, types, has_dorsal, has_ventral)
% Returns a struct array (one entry per trial type) with N, mean, SEM of
% the per-session-normalized SWR rates for dorsal and ventral probes.
% Total counts are omitted (they are in computeSummary).

summary = struct();
for i = 1:numel(types)
    summary(i).type = types{i};
    idx = strcmp({ITI_results.type_before}, types{i});

    if has_dorsal
        vals = [ITI_results(idx).rate_dorsal_SWR_norm];
        n    = sum(~isnan(vals));
        summary(i).dorsal_N    = n;
        summary(i).dorsal_mean = mean(vals, 'omitnan');
        summary(i).dorsal_sem  = std(vals, 'omitnan') / max(sqrt(n), 1);
    else
        summary(i).dorsal_N    = NaN;
        summary(i).dorsal_mean = NaN;
        summary(i).dorsal_sem  = NaN;
    end

    if has_ventral
        vals = [ITI_results(idx).rate_ventral_SWR_norm];
        n    = sum(~isnan(vals));
        summary(i).ventral_N    = n;
        summary(i).ventral_mean = mean(vals, 'omitnan');
        summary(i).ventral_sem  = std(vals, 'omitnan') / max(sqrt(n), 1);
    else
        summary(i).ventral_N    = NaN;
        summary(i).ventral_mean = NaN;
        summary(i).ventral_sem  = NaN;
    end
end
end

% -------------------------------------------------------------------------
function summary = computeSummaryNormHierarchical(all_sessions, types, has_dorsal, has_ventral)
% Computes combined normalized SWR rate using two-level (hierarchical) averaging.
%
% Level 1: each session already has a per-trial-type mean in summary_norm
%           (computed over that session's ITIs).
% Level 2: then we average those session-level means across sessions,
%           so that every session contributes equally regardless of ITI count.
%
% The SEM thus reflects session-to-session variabilit

    n_sessions = numel(all_sessions);
    summary    = struct();
    
    for i = 1:numel(types)
        summary(i).type = types{i};
    
        sess_means_dorsal  = NaN(n_sessions, 1);
        sess_means_ventral = NaN(n_sessions, 1);
    
        for s = 1:n_sessions
            sess_norm = all_sessions(s).summary_norm;          % per-session summary_norm struct
    
            % find the matching trial type entry in this session's summary
            type_match = strcmp({sess_norm.type}, types{i});
            if ~any(type_match)
                continue                                        % trial type absent in this session
            end
    
            if has_dorsal && all_sessions(s).has_dorsal
                n_iterable = sess_norm(type_match).dorsal_N;
                if ~isnan(n_iterable) && n_iterable > 0
                    sess_means_dorsal(s) = sess_norm(type_match).dorsal_mean;
                end
            end
    
            if has_ventral && all_sessions(s).has_ventral
                n_iterable = sess_norm(type_match).ventral_N;
                if ~isnan(n_iterable) && n_iterable > 0
                    sess_means_ventral(s) = sess_norm(type_match).ventral_mean;
                end
            end
        end
    
        % --- Dorsal ---
        if has_dorsal
            valid = sess_means_dorsal(~isnan(sess_means_dorsal));
            n     = numel(valid);
            summary(i).dorsal_N    = n;
            summary(i).dorsal_mean = mean(valid);
            summary(i).dorsal_sem  = std(valid) / max(sqrt(n), 1);
            summary(i).dorsal_session_means = sess_means_dorsal; % keep for inspection / stats
        else
            summary(i).dorsal_N    = NaN;
            summary(i).dorsal_mean = NaN;
            summary(i).dorsal_sem  = NaN;
            summary(i).dorsal_session_means = NaN;
        end
    
        % --- Ventral ---
        if has_ventral
            valid = sess_means_ventral(~isnan(sess_means_ventral));
            n     = numel(valid);
            summary(i).ventral_N    = n;
            summary(i).ventral_mean = mean(valid);
            summary(i).ventral_sem  = std(valid) / max(sqrt(n), 1);
            summary(i).ventral_session_means = sess_means_ventral;
        else
            summary(i).ventral_N    = NaN;
            summary(i).ventral_mean = NaN;
            summary(i).ventral_sem  = NaN;
            summary(i).ventral_session_means = NaN;
        end
    end
end

% -------------------------------------------------------------------------
function fig = plotSWRbarsNormHierarchical(summary, types, has_dorsal, has_ventral, fig_title)
% Bar plot for the hierarchical normalized summary.
% Each bar = mean across sessions; error bar = SEM across sessions.
% Individual session means are overlaid as dots so you can see the spread.

    n_probes = has_dorsal + has_ventral;
    fig = figure('Name', fig_title, 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 max(n_probes*14, 14) 12]);
    colors  = lines(numel(types));
    sp_idx  = 0;
    
    probe_flags    = [has_dorsal,       has_ventral];
    mean_fields    = {'dorsal_mean',    'ventral_mean'};
    sem_fields     = {'dorsal_sem',     'ventral_sem'};
    n_fields       = {'dorsal_N',       'ventral_N'};
    sess_fields    = {'dorsal_session_means', 'ventral_session_means'};
    titles         = {'Dorsal normalized SWR rate in ITI', ...
                      'Ventral normalized SWR rate in ITI'};
    
    for p = 1:2
        if ~probe_flags(p), continue; end
        sp_idx = sp_idx + 1;
        ax = subplot(1, n_probes, sp_idx);
        hold(ax, 'on');
    
        means = zeros(1, numel(types));
        sems  = zeros(1, numel(types));
        ns    = zeros(1, numel(types));
    
        for i = 1:numel(types)
            means(i) = summary(i).(mean_fields{p});
            sems(i)  = summary(i).(sem_fields{p});
            ns(i)    = summary(i).(n_fields{p});
        end
    
        % Replace NaN with 0 just for the bar call so MATLAB doesn't skip bars
        means_plot = means;
        means_plot(isnan(means_plot)) = 0;
    
        bar(ax, 1:numel(types), means_plot, 'FaceColor', 'flat', 'CData', colors);
        errorbar(ax, 1:numel(types), means, sems, 'k.', 'LineWidth', 1);
    
        % Overlay individual session means as dots
        for i = 1:numel(types)
            sess_vals = summary(i).(sess_fields{p});
            sess_vals = sess_vals(~isnan(sess_vals));
            if isempty(sess_vals), continue; end
            jitter = (rand(size(sess_vals)) - 0.5) * 0.2;   % small horizontal jitter
            scatter(ax, i + jitter, sess_vals, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
        end
    
        set(ax, 'XTick', 1:numel(types), 'XTickLabel', types, ...
            'XTickLabelRotation', 30, 'FontSize', 11);
        ylabel(ax, 'Normalized SWR rate (% of session / s)', 'FontSize', 14);
        title(ax, titles{p}, 'FontSize', 15);
        box(ax, 'off');
    
        % Label number of sessions above each bar
        y_max    = max(means + sems, [], 'omitnan');
        y_offset = y_max * 0.01;
        for i = 1:numel(types)
            if ns(i) > 0
                text(ax, i, means(i) + sems(i) + y_offset, ...
                    sprintf('n sess=%d', ns(i)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'bottom', ...
                    'FontSize', 6);
            end
        end
    end
    
    sgtitle(fig_title, 'FontSize', 14, 'Interpreter', 'none');
end

% -------------------------------------------------------------------------
function printSummary(summary, has_dorsal, has_ventral)
% Prints a formatted summary table (raw rates) to the command window.

    if has_dorsal
        fprintf('\n  Dorsal SWR rate (SWR/s) by preceding trial type:\n');
        fprintf('  %-25s  %6s  %8s  %8s  %10s\n', ...
            'Trial type (before)', 'N ITIs', 'Mean', 'SEM', 'Total SWRs');
        for i = 1:numel(summary)
            if isnan(summary(i).dorsal_N) || summary(i).dorsal_N == 0, continue; end
            fprintf('  %-25s  %6d  %8.4f  %8.4f  %10d\n', ...
                summary(i).type, summary(i).dorsal_N, ...
                summary(i).dorsal_mean, summary(i).dorsal_sem, ...
                summary(i).dorsal_total);
        end
    end
    
    if has_ventral
        fprintf('\n  Ventral SWR rate (SWR/s) by preceding trial type:\n');
        fprintf('  %-25s  %6s  %8s  %8s  %10s\n', ...
            'Trial type (before)', 'N ITIs', 'Mean', 'SEM', 'Total SWRs');
        for i = 1:numel(summary)
            if isnan(summary(i).ventral_N) || summary(i).ventral_N == 0, continue; end
            fprintf('  %-25s  %6d  %8.4f  %8.4f  %10d\n', ...
                summary(i).type, summary(i).ventral_N, ...
                summary(i).ventral_mean, summary(i).ventral_sem, ...
                summary(i).ventral_total);
        end
    end
end

% -------------------------------------------------------------------------
function printSummaryNorm(summary, has_dorsal, has_ventral)
% Prints a formatted summary table (normalized rates) to the command window.

    if has_dorsal
        fprintf('\n  Dorsal normalized SWR rate (percentage of session SWRs / s) by preceding trial type:\n');
        fprintf('  %-25s  %6s  %8s  %8s\n', ...
            'Trial type (before)', 'N ITIs', 'Mean', 'SEM');
        for i = 1:numel(summary)
            if isnan(summary(i).dorsal_N) || summary(i).dorsal_N == 0, continue; end
            fprintf('  %-25s  %6d  %8.6f  %8.6f\n', ...
                summary(i).type, summary(i).dorsal_N, ...
                summary(i).dorsal_mean, summary(i).dorsal_sem);
        end
    end
    
    if has_ventral
        fprintf('\n  Ventral normalized SWR rate (percentage of session SWRs / s) by preceding trial type:\n');
        fprintf('  %-25s  %6s  %8s  %8s\n', ...
            'Trial type (before)', 'N ITIs', 'Mean', 'SEM');
        for i = 1:numel(summary)
            if isnan(summary(i).ventral_N) || summary(i).ventral_N == 0, continue; end
            fprintf('  %-25s  %6d  %8.6f  %8.6f\n', ...
                summary(i).type, summary(i).ventral_N, ...
                summary(i).ventral_mean, summary(i).ventral_sem);
        end
    end
end

% -------------------------------------------------------------------------
function fig = plotSWRbars(ITI_results, types, has_dorsal, has_ventral, fig_title)
% Generates a bar figure of raw SWR rates (one subplot per probe).
% Labels above bars show total raw SWR count.

    n_probes = has_dorsal + has_ventral;
    fig = figure('Name', fig_title, 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 max(n_probes*14, 14) 12]);
    colors  = lines(numel(types));
    sp_idx  = 0;
    
    probe_flags = [has_dorsal, has_ventral];
    rate_fields = {'rate_dorsal_SWR',  'rate_ventral_SWR'};
    n_fields    = {'n_dorsal_SWR',     'n_ventral_SWR'};
    titles      = {'Dorsal SWR rate in ITI', 'Ventral SWR rate in ITI'};
    
    for p = 1:2
        if ~probe_flags(p), continue; end
        sp_idx = sp_idx + 1;
        ax = subplot(1, n_probes, sp_idx);
        hold(ax, 'on');
    
        means  = zeros(1, numel(types));
        sems   = zeros(1, numel(types));
        totals = zeros(1, numel(types));
    
        for i = 1:numel(types)
            idx  = strcmp({ITI_results.type_before}, types{i});
            vals = [ITI_results(idx).(rate_fields{p})];
            if isempty(vals) || all(isnan(vals)), continue; end
            n         = sum(~isnan(vals));
            means(i)  = mean(vals, 'omitnan');
            sems(i)   = std(vals,  'omitnan') / max(sqrt(n), 1);
            totals(i) = sum([ITI_results(idx).(n_fields{p})], 'omitnan');
        end
    
        bar(ax, 1:numel(types), means, 'FaceColor', 'flat', 'CData', colors);
        errorbar(ax, 1:numel(types), means, sems, 'k.', 'LineWidth', 1);
        set(ax, 'XTick', 1:numel(types), 'XTickLabel', types, ...
            'XTickLabelRotation', 30, 'FontSize', 11);
        ylabel(ax, 'SWR rate (events/s)', 'FontSize', 14);
        title(ax, titles{p}, 'FontSize', 15);
        box(ax, 'off');
    
        y_max    = max(means + sems);
        y_offset = y_max * 0.01;
        for i = 1:numel(types)
            if totals(i) > 0
                text(ax, i, means(i) + sems(i) + y_offset, ...
                    sprintf('n SWRs=%d', totals(i)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'bottom', ...
                    'FontSize', 6);
            end
        end
    end
    
    sgtitle(fig_title, 'FontSize', 14, 'Interpreter', 'none');
end

% -------------------------------------------------------------------------
function fig = plotSWRbarsNorm(ITI_results, types, has_dorsal, has_ventral, fig_title)
% Generates a bar figure of normalized SWR rates (one subplot per probe).
% Each ITI's rate = n_SWR / (session_total_SWR * duration), so bars
% represent each trial type's share of the session's SWR per second.
% Labels above bars show the number of ITIs contributing to each bar.

    n_probes = has_dorsal + has_ventral;
    fig = figure('Name', fig_title, 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 max(n_probes*14, 14) 12]);
    colors  = lines(numel(types));
    sp_idx  = 0;
    
    probe_flags  = [has_dorsal, has_ventral];
    rate_fields  = {'rate_dorsal_SWR_norm', 'rate_ventral_SWR_norm'};
    titles       = {'Dorsal normalized SWR rate in ITI', 'Ventral normalized SWR rate in ITI'};
    
    for p = 1:2
        if ~probe_flags(p), continue; end
        sp_idx = sp_idx + 1;
        ax = subplot(1, n_probes, sp_idx);
        hold(ax, 'on');
    
        means = zeros(1, numel(types));
        sems  = zeros(1, numel(types));
        ns    = zeros(1, numel(types));
    
        for i = 1:numel(types)
            idx  = strcmp({ITI_results.type_before}, types{i});
            vals = [ITI_results(idx).(rate_fields{p})];
            if isempty(vals) || all(isnan(vals)), continue; end
            n        = sum(~isnan(vals));
            means(i) = mean(vals, 'omitnan');
            sems(i)  = std(vals,  'omitnan') / max(sqrt(n), 1);
            ns(i)    = n;
        end
    
        bar(ax, 1:numel(types), means, 'FaceColor', 'flat', 'CData', colors);
        errorbar(ax, 1:numel(types), means, sems, 'k.', 'LineWidth', 1);
        set(ax, 'XTick', 1:numel(types), 'XTickLabel', types, ...
            'XTickLabelRotation', 30, 'FontSize', 11);
        ylabel(ax, 'Normalized SWR rate (% of session / s)', 'FontSize', 14);
        title(ax, titles{p}, 'FontSize', 15);
        box(ax, 'off');
    
        % Label: number of ITIs above each bar
        y_max    = max(means + sems);
        y_offset = y_max * 0.01;
        for i = 1:numel(types)
            if ns(i) > 0
                text(ax, i, means(i) + sems(i) + y_offset, ...
                    sprintf('n ITIs=%d', ns(i)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'bottom', ...
                    'FontSize', 6);
            end
        end
    end
    
    sgtitle(fig_title, 'FontSize', 14, 'Interpreter', 'none');
end

% -------------------------------------------------------------------------
function [tstart, tend] = extractSWRtimes(s, struct_name)
% Extracts SWR start/end times from s.iv.tstart and s.iv.tend.

    if ~isfield(s, 'iv')
        error('Expected field ".iv" not found in %s', struct_name);
    end
    iv = s.iv;
    if ~isfield(iv, 'tstart') || ~isfield(iv, 'tend')
        error('Expected ".iv.tstart" / ".iv.tend" not found in %s', struct_name);
    end
    tstart = double(iv.tstart(:));
    tend   = double(iv.tend(:));
    if numel(tstart) ~= numel(tend)
        error('tstart and tend have different lengths in %s', struct_name);
    end
    fprintf('Loaded %d SWRs from "%s"\n', numel(tstart), struct_name);
end

% -------------------------------------------------------------------------
function times = getEventTimes(ev, label)
% Returns event onset times matching a label string (case-insensitive).

    times = [];
    if ~(isfield(ev, 'label') && isfield(ev, 't'))
        warning('no event times for "%s" found in all_events', label);
        return;
    end
    labels = ev.label;
    for k = 1:numel(labels)
        if strcmpi(strtrim(labels{k}), strtrim(label))
            t = ev.t;
            if iscell(t)
                times = double(t{k}(:));
            else
                times = double(t(k,:));
                times = times(~isnan(times));
            end
            return;
        end
    end
end

% -------------------------------------------------------------------------
function trials = classifyTrials( ...
    rew_od1, ap_od1, neu_od1, ...
    rew_out, ap_out, ...
    win_min, win_max, catch_win, neutral_win, rew_trial_buff, ap_trial_buff)
    % Classify each event occurrence into a named trial type.
    
    trials = struct('type', {}, 'startTime', {}, 'endTime', {});
    idx = 0;
    
    % ---- 1. REWARD trials ---------------------------------------------------
    for k = 1:numel(rew_od1)
        t0    = rew_od1(k);
        out_t = rew_out(rew_out >= t0+win_min & rew_out <= t0+win_max);
        idx   = idx + 1;
        if ~isempty(out_t)
            trials(idx).type      = 'reward';
            trials(idx).startTime = t0;
            trials(idx).endTime   = out_t(1) + rew_trial_buff;
        else
            trials(idx).type      = 'odor-reward-catch';
            trials(idx).startTime = t0;
            trials(idx).endTime   = t0 + catch_win;
        end
    end
    
    % ---- 2. AIRPUFF trials --------------------------------------------------
    for k = 1:numel(ap_od1)
        t0    = ap_od1(k);
        out_t = ap_out(ap_out >= t0+win_min & ap_out <= t0+win_max);
        idx   = idx + 1;
        if ~isempty(out_t)
            trials(idx).type      = 'airpuff';
            trials(idx).startTime = t0;
            trials(idx).endTime   = out_t(1) + ap_trial_buff;
        else
            trials(idx).type      = 'odor-airpuff-catch';
            trials(idx).startTime = t0;
            trials(idx).endTime   = t0 + catch_win;
        end
    end
    
    % ---- 3. NEUTRAL trials --------------------------------------------------
    for k = 1:numel(neu_od1)
        t0  = neu_od1(k);
        idx = idx + 1;
        trials(idx).type      = 'neutral';
        trials(idx).startTime = t0;
        trials(idx).endTime   = t0 + neutral_win;
    end
    
    % ---- 4. REWARD-CATCH ----------------------------------------------------
    for k = 1:numel(rew_out)
        t_out    = rew_out(k);
        preceded = any(rew_od1 >= t_out-win_max & rew_od1 <= t_out-win_min);
        if ~preceded
            idx = idx + 1;
            trials(idx).type      = 'reward-catch';
            trials(idx).startTime = t_out;
            trials(idx).endTime   = t_out + rew_trial_buff;
        end
    end
    
    % ---- 5. AIRPUFF-CATCH ---------------------------------------------------
    for k = 1:numel(ap_out)
        t_out    = ap_out(k);
        preceded = any(ap_od1 >= t_out-win_max & ap_od1 <= t_out-win_min);
        if ~preceded
            idx = idx + 1;
            trials(idx).type      = 'airpuff-catch';
            trials(idx).startTime = t_out;
            trials(idx).endTime   = t_out + ap_trial_buff;
        end
    end
end

