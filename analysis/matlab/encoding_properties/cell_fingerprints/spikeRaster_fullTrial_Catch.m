%% SPIKE RASTER PLOT GENERATOR - FULL TRIAL VERSION (CATCH TRIAL ADAPTED)
% This script generates raster plots for neural spiking data across full trials
% Each plot shows spike timing for a neuron across multiple trials of the same type
% (Reward trials, Airpuff trials, or Neutral trials)
%
% CATCH TRIAL HANDLING:
%   - Normal trials (odor pair + US): black spikes, top of raster
%   - Odor-only catch (odor pair, no US): orange spikes, middle of raster
%   - US-only catch (US, no odor pair): purple spikes, aligned to US, bottom of raster
%   Sections are separated by gaps and labeled on the y-axis
%
% TRIAL TIMING (from odor 1 onset):
%   - Odor 1 onset:        0.00 s
%   - Odor 1 offset:       1.00 s
%   - Break:               1.00 - 1.50 s
%   - Odor 2 onset:        1.50 s
%   - Odor 2 offset:       2.50 s
%   - Airpuff delivery:    2.75 s (0.25 s after odor 2 offset)
%   - Reward delivery:     3.00 s (0.50 s after odor 2 offset)

%% Input information
cd('D:\vCA1\M649\preprocessed\M649_2026_04_01');
cfg.probe = 1;
cfg.spikeIDs = {'imec0_341'};
cfg.pre_trial = 1.0;   % time before odor 1 onset to plot (seconds)
cfg.post_trial = 5.0;  % time after outcome to plot (seconds)

addpath(genpath('C:\Users\mvdmlab\Documents\egibbens\odor-pixels'));
LoadExpKeys;

%% Load Spikes
if cfg.probe == 1
    S = load_NPX_spike('clean_units_imec0.mat');
    record_depth = ExpKeys.probe1_Depth;
    cfg_md.maxdepth = 1000 * record_depth;
    S.usr.recording_depths = (cfg_md.maxdepth - S.usr.recording_depths) * cos(deg2rad(ExpKeys.probe1_insertion_roll));
else
    S = load_NPX_spike('clean_units_imec1.mat');
    record_depth = ExpKeys.probe2_Depth;
    cfg_md.maxdepth = 1000 * record_depth;
    S.usr.recording_depths = (cfg_md.maxdepth - S.usr.recording_depths) * cos(deg2rad(ExpKeys.probe2_insertion_roll));
end

% Sort all units according to depth
[~, sort_idx] = sort(S.usr.recording_depths, 'descend');
S.t = S.t(sort_idx);
S.label = S.label(sort_idx);
S.usr.mean_waveforms = S.usr.mean_waveforms(sort_idx);
S.usr.recording_depths = S.usr.recording_depths(sort_idx);
S.usr.shank_ids = S.usr.shank_ids(sort_idx);
clear sort_idx;

%% Load events and extract stimulus times
load('all_events.mat')

% Extract all odor onsets and identities together
all_odor_onsets = [];
all_odor_identities = {};

% Get all odor names from ExpKeys
all_possible_odors = [ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors];

% Extract each odor's onsets
for odor_idx = 1:length(all_possible_odors)
    odor_name = all_possible_odors{odor_idx};
    odor_mask = strcmp(evt.label, odor_name);

    if any(odor_mask)
        odor_times = evt.t(odor_mask);
        if iscell(odor_times) && ~isempty(odor_times{1})
            times = odor_times{1};
            if ~isempty(times)
                all_odor_onsets = [all_odor_onsets; times(:)];
                all_odor_identities = [all_odor_identities; repmat({odor_name}, length(times), 1)];
            end
        end
    end
end

% Sort by time
if ~isempty(all_odor_onsets)
    [odor_onsets, sort_idx] = sort(all_odor_onsets);
    odor_identities = all_odor_identities(sort_idx);
else
    error('No odor events found in the event file.');
end

% Extract unconditioned stimulus times
reward_times = evt.t(strcmp(evt.label, 'Reward ON'));
if ~isempty(reward_times)
    if iscell(reward_times)
        reward_times = [reward_times{:}];
    end
else
    reward_times = [];
end

airpuff_times = evt.t(strcmp(evt.label, 'Airpuff ON'));
if ~isempty(airpuff_times)
    if iscell(airpuff_times)
        airpuff_times = [airpuff_times{:}];
    end
else
    airpuff_times = [];
end

%% Classify odors into categories using ExpKeys
odor_to_outcome = containers.Map();

if isfield(ExpKeys, 'reward_odors')
    for i = 1:length(ExpKeys.reward_odors)
        odor_to_outcome(ExpKeys.reward_odors{i}) = 'Reward';
    end
end

if isfield(ExpKeys, 'airpuff_odors')
    for i = 1:length(ExpKeys.airpuff_odors)
        odor_to_outcome(ExpKeys.airpuff_odors{i}) = 'Airpuff';
    end
end

if isfield(ExpKeys, 'neutral_odors')
    for i = 1:length(ExpKeys.neutral_odors)
        odor_to_outcome(ExpKeys.neutral_odors{i}) = 'Neutral';
    end
end

%% Fixed US latencies from odor 1 onset
airpuff_latency_from_odor1 = 2.75;  % seconds
reward_latency_from_odor1  = 3.00;  % seconds

%% Identify odor 1 labels for each trial category
% First odor in each ExpKeys pair is odor 1
reward_odor1_label  = ExpKeys.reward_odors{1};
reward_odor2_label  = ExpKeys.reward_odors{2};
airpuff_odor1_label = ExpKeys.airpuff_odors{1};
airpuff_odor2_label = ExpKeys.airpuff_odors{2};
neutral_odor1_label = ExpKeys.neutral_odors{1};
neutral_odor2_label = ExpKeys.neutral_odors{2};

%% Build trial structure — classify normal vs catch

% Define odor pairs: {odor1_label, odor2_label, outcome_type, expected_US_latency}
pair_defs = {
    reward_odor1_label,  reward_odor2_label,  'Reward',  reward_latency_from_odor1;
    airpuff_odor1_label, airpuff_odor2_label, 'Airpuff', airpuff_latency_from_odor1;
    neutral_odor1_label, neutral_odor2_label, 'Neutral', NaN;
    };

pair_tol = 0.5;  % tolerance for matching US to odor pair (seconds)
max_odor_pair_gap = 5;  % max seconds between odor 1 and odor 2

trials = struct('odor1_onset', {}, 'odor1_name', {}, 'odor2_onset', {}, ...
    'odor2_name', {}, 'outcome_time', {}, 'outcome_type', {}, 'catch_type', {});
trial_count = 0;

% Track which US events get claimed by a paired trial
reward_claimed  = false(size(reward_times));
airpuff_claimed = false(size(airpuff_times));

for pair_idx = 1:size(pair_defs, 1)
    odor1_name    = pair_defs{pair_idx, 1};
    odor2_name    = pair_defs{pair_idx, 2};
    outcome_type  = pair_defs{pair_idx, 3};
    expected_lat  = pair_defs{pair_idx, 4};

    % Find all onsets of odor 1 for this pair
    odor1_mask = strcmp(odor_identities, odor1_name);
    odor1_onsets_this = odor_onsets(odor1_mask);

    for i = 1:length(odor1_onsets_this)
        odor1_time = odor1_onsets_this(i);

        % Find the next odor 2 onset after this odor 1
        subsequent_odor2 = odor_onsets(strcmp(odor_identities, odor2_name) & odor_onsets > odor1_time);

        if ~isempty(subsequent_odor2)
            odor2_time = subsequent_odor2(1);

            % Check valid pair (odor 2 close enough in time)
            if odor2_time - odor1_time < max_odor_pair_gap
                trial_count = trial_count + 1;
                trials(trial_count).odor1_onset = odor1_time;
                trials(trial_count).odor1_name  = odor1_name;
                trials(trial_count).odor2_onset = odor2_time;
                trials(trial_count).odor2_name  = odor2_name;
                trials(trial_count).outcome_type = outcome_type;

                % Try to find a matching US
                switch outcome_type
                    case 'Reward'
                        dt = reward_times - odor1_time;               % time from odor‑1 to each reward
                        % make sure the claim‑vector has the same orientation as dt
                        %   (reward_claimed is already a column vector of length reward_times)
                        valid = dt > 0 & abs(dt - expected_lat) <= pair_tol & ~reward_claimed;

                        if any(valid)                     % <‑‑ protect against an empty set
                            % find the reward that is closest to the expected latency
                            [~, best] = min(abs(dt(valid) - expected_lat));
                            valid_idx = find(valid);      % indices of all “eligible” rewards
                            us_idx = valid_idx(best);     % the one we will claim
                            trials(trial_count).outcome_time = reward_times(us_idx);
                            trials(trial_count).catch_type   = 'normal';
                            reward_claimed(us_idx) = true;   % mark this reward as used
                        else
                            % No matching reward → this is an odor‑only catch trial
                            trials(trial_count).outcome_time = odor1_time + expected_lat;
                            trials(trial_count).catch_type   = 'odor_only_catch';
                        end

                    case 'Airpuff'
                        dt = airpuff_times - odor1_time;
                        valid = dt > 0 & abs(dt - expected_lat) <= pair_tol & ~airpuff_claimed;

                        if any(valid)
                            [~, best] = min(abs(dt(valid) - expected_lat));
                            valid_idx = find(valid);
                            us_idx = valid_idx(best);
                            trials(trial_count).outcome_time = airpuff_times(us_idx);
                            trials(trial_count).catch_type   = 'normal';
                            airpuff_claimed(us_idx) = true;
                        else
                            trials(trial_count).outcome_time = odor1_time + expected_lat;
                            trials(trial_count).catch_type   = 'odor_only_catch';
                        end

                    case 'Neutral'
                        trials(trial_count).outcome_time = odor1_time + 3.5;  % no US expected
                        trials(trial_count).catch_type = 'normal';
                end
            end
        end
    end
end

% Build US-only catch trials from unclaimed US events
% Reward US-only catches
unclaimed_reward = reward_times(~reward_claimed);
for i = 1:length(unclaimed_reward)
    trial_count = trial_count + 1;
    trials(trial_count).odor1_onset  = NaN;
    trials(trial_count).odor1_name   = '';
    trials(trial_count).odor2_onset  = NaN;
    trials(trial_count).odor2_name   = '';
    trials(trial_count).outcome_time = unclaimed_reward(i);
    trials(trial_count).outcome_type = 'Reward';
    trials(trial_count).catch_type   = 'us_only_catch';
end

% Airpuff US-only catches
unclaimed_airpuff = airpuff_times(~airpuff_claimed);
for i = 1:length(unclaimed_airpuff)
    trial_count = trial_count + 1;
    trials(trial_count).odor1_onset  = NaN;
    trials(trial_count).odor1_name   = '';
    trials(trial_count).odor2_onset  = NaN;
    trials(trial_count).odor2_name   = '';
    trials(trial_count).outcome_time = unclaimed_airpuff(i);
    trials(trial_count).outcome_type = 'Airpuff';
    trials(trial_count).catch_type   = 'us_only_catch';
end

% Print diagnostics
fprintf('\n=== Trial classification ===\n');
for type = {'Reward', 'Airpuff', 'Neutral'}
    n_normal    = sum(strcmp({trials.outcome_type}, type{1}) & strcmp({trials.catch_type}, 'normal'));
    n_odor_catch = sum(strcmp({trials.outcome_type}, type{1}) & strcmp({trials.catch_type}, 'odor_only_catch'));
    n_us_catch   = sum(strcmp({trials.outcome_type}, type{1}) & strcmp({trials.catch_type}, 'us_only_catch'));
    fprintf('%s: %d normal, %d odor-only catch, %d US-only catch\n', ...
        type{1}, n_normal, n_odor_catch, n_us_catch);
end
fprintf('============================\n\n');

%% Isolate spike data for neurons of interest
S_plot = S;
if ~isempty(cfg.spikeIDs)
    spikeIdx = find(cellfun(@(x) ismember(x, cfg.spikeIDs), S.label));
    if isempty(spikeIdx)
        error('No neurons found matching cfg.spikeIDs. Check neuron labels.');
    end
    S_plot.label = S.label(spikeIdx);
    S_plot.t = S.t(spikeIdx);
    S_plot.usr.channel_ids = S.usr.channel_ids(spikeIdx);
    S_plot.usr.recording_depths = S.usr.recording_depths(spikeIdx);
    S_plot.usr.shank_ids = S.usr.shank_ids(spikeIdx);
    S_plot.usr.mean_waveforms = S.usr.mean_waveforms(spikeIdx);
end

%% Generate full trial raster plots for each neuron
if ~isempty(trials)
    generate_full_trial_rasters_catch(S_plot, trials, cfg.pre_trial, cfg.post_trial, ...
        reward_latency_from_odor1, airpuff_latency_from_odor1, './full_trial_rasters');
else
    warning('No valid trials found. Check event file and odor labels.');
end

%% Function to generate full trial raster plots with catch trial support
function generate_full_trial_rasters_catch(S, trials, time_before, time_after, ...
    reward_lat, airpuff_lat, output_dir)
% GENERATE_FULL_TRIAL_RASTERS_CATCH Creates raster plots with catch trial support
%
% Trial sections (top to bottom, separated by gaps):
%   1. Normal trials (black spikes)
%   2. Odor-only catch trials (orange spikes) - aligned to odor 1 onset
%   3. US-only catch trials (purple spikes) - aligned to US onset
%
% Normal and odor-only catch trials are aligned to odor 1 onset (t=0).
% US-only catch trials are aligned to US onset (t=0) since there is no odor.

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Define trial types for separate figures
trial_types = {'Reward', 'Airpuff', 'Neutral'};

% Spike colors for each section
color_normal     = [0 0 0];          % black
color_odor_catch = [0.85 0.45 0];    % orange
color_us_catch   = [0.5 0 0.8];     % purple

% Colors for outcome reference lines
outcome_colors = containers.Map({'Reward', 'Airpuff', 'Neutral'}, ...
    {[0, 0, 1], [1, 0, 0], [0.5, 0.5, 0.5]});

% Gap between sections (in trial units)
section_gap = 2;

% Loop through each neuron
for neuron_idx = 1:length(S.t)
    spike_times = S.t{neuron_idx};
    neuron_label = S.label{neuron_idx};

    % Create one figure per trial type
    for trial_type_idx = 1:length(trial_types)
        current_trial_type = trial_types{trial_type_idx};
        outcome_color = outcome_colors(current_trial_type);

        % Get expected US latency for this trial type
        switch current_trial_type
            case 'Reward'
                expected_us_lat = reward_lat;
            case 'Airpuff'
                expected_us_lat = airpuff_lat;
            case 'Neutral'
                expected_us_lat = NaN;
        end

        % Separate trials into sections
        idx_normal     = find(strcmp({trials.outcome_type}, current_trial_type) & ...
            strcmp({trials.catch_type}, 'normal'));
        idx_odor_catch = find(strcmp({trials.outcome_type}, current_trial_type) & ...
            strcmp({trials.catch_type}, 'odor_only_catch'));
        idx_us_catch   = find(strcmp({trials.outcome_type}, current_trial_type) & ...
            strcmp({trials.catch_type}, 'us_only_catch'));

        n_normal     = length(idx_normal);
        n_odor_catch = length(idx_odor_catch);
        n_us_catch   = length(idx_us_catch);

        % Skip if no trials at all for this type
        if n_normal == 0 && n_odor_catch == 0 && n_us_catch == 0
            fprintf('No %s trials found for %s\n', current_trial_type, neuron_label);
            continue;
        end

        % Calculate y-axis positions for each section
        % Section 1: Normal trials (top)
        normal_y_start = 1;
        normal_y_end   = n_normal;

        % Section 2: Odor-only catch (middle)
        if n_odor_catch > 0
            odor_catch_y_start = normal_y_end + section_gap + 1;
        else
            odor_catch_y_start = normal_y_end + 1;
        end
        odor_catch_y_end = odor_catch_y_start + n_odor_catch - 1;

        % Section 3: US-only catch (bottom)
        if n_us_catch > 0
            if n_odor_catch > 0
                us_catch_y_start = odor_catch_y_end + section_gap + 1;
            else
                us_catch_y_start = normal_y_end + section_gap + 1;
            end
        else
            us_catch_y_start = odor_catch_y_end + 1;
        end
        us_catch_y_end = us_catch_y_start + n_us_catch - 1;

        % Total height of the raster
        max_y = max([normal_y_end, odor_catch_y_end, us_catch_y_end, 1]);

        % Create figure
        fig = figure('Position', [100, 100, 1600, 900]);
        hold on;

        % =============================================================
        % SECTION 1: Normal trials (aligned to odor 1 onset, black spikes)
        % =============================================================
        for plot_idx = 1:n_normal
            trial_num = idx_normal(plot_idx);
            trial = trials(trial_num);
            y_pos = normal_y_start + plot_idx - 1;

            % Define window around odor 1 onset
            trial_start = trial.odor1_onset - time_before;
            trial_end   = trial.odor1_onset + max([reward_lat, airpuff_lat]) + time_after;

            % Find spikes within window
            spikes_in_trial = spike_times(spike_times >= trial_start & spike_times <= trial_end);
            relative_spikes = spikes_in_trial - trial.odor1_onset;

            % Plot spikes (black)
            for spike = relative_spikes'
                plot([spike, spike], [y_pos - 0.4, y_pos + 0.4], '-', ...
                    'Color', color_normal, 'LineWidth', 0.5);
            end

            % Per-trial event markers
            % Odor 1 onset (red solid)
            line([0, 0], [y_pos - 0.4, y_pos + 0.4], ...
                'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);

            % Odor 2 onset (blue dotted)
            odor2_relative = trial.odor2_onset - trial.odor1_onset;
            line([odor2_relative, odor2_relative], [y_pos - 0.4, y_pos + 0.4], ...
                'Color', 'b', 'LineStyle', ':', 'LineWidth', 1.5);

            % Outcome marker (dashed, outcome color)
            outcome_relative = trial.outcome_time - trial.odor1_onset;
            line([outcome_relative, outcome_relative], [y_pos - 0.4, y_pos + 0.4], ...
                'Color', outcome_color, 'LineStyle', '--', 'LineWidth', 1.5);
        end

        % =============================================================
        % SECTION 2: Odor-only catch trials (aligned to odor 1 onset, orange spikes)
        % =============================================================
        for plot_idx = 1:n_odor_catch
            trial_num = idx_odor_catch(plot_idx);
            trial = trials(trial_num);
            y_pos = odor_catch_y_start + plot_idx - 1;

            % Define window around odor 1 onset
            trial_start = trial.odor1_onset - time_before;
            trial_end   = trial.odor1_onset + max([reward_lat, airpuff_lat]) + time_after;

            % Find spikes within window
            spikes_in_trial = spike_times(spike_times >= trial_start & spike_times <= trial_end);
            relative_spikes = spikes_in_trial - trial.odor1_onset;

            % Plot spikes (orange)
            for spike = relative_spikes'
                plot([spike, spike], [y_pos - 0.4, y_pos + 0.4], '-', ...
                    'Color', color_odor_catch, 'LineWidth', 0.5);
            end

            % Odor 1 onset (red solid)
            line([0, 0], [y_pos - 0.4, y_pos + 0.4], ...
                'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);

            % Odor 2 onset (blue dotted)
            odor2_relative = trial.odor2_onset - trial.odor1_onset;
            line([odor2_relative, odor2_relative], [y_pos - 0.4, y_pos + 0.4], ...
                'Color', 'b', 'LineStyle', ':', 'LineWidth', 1.5);

            % Expected outcome time (dashed gray — no actual US delivered)
            if ~isnan(expected_us_lat)
                line([expected_us_lat, expected_us_lat], [y_pos - 0.4, y_pos + 0.4], ...
                    'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.5);
            end
        end

        % =============================================================
        % SECTION 3: US-only catch trials (aligned to US onset, purple spikes)
        % =============================================================
        for plot_idx = 1:n_us_catch
            trial_num = idx_us_catch(plot_idx);
            trial = trials(trial_num);
            y_pos = us_catch_y_start + plot_idx - 1;

            % For US-only catch, align to US onset
            us_time = trial.outcome_time;
            trial_start = us_time - time_before - max([reward_lat, airpuff_lat]);
            trial_end   = us_time + time_after;

            % Find spikes within window
            spikes_in_trial = spike_times(spike_times >= trial_start & spike_times <= trial_end);

            % Align spikes to where US would be in the normal trial timeline
            % Place US at the expected latency position so it lines up with
            % the normal trial outcome markers above
            relative_spikes = spikes_in_trial - us_time + expected_us_lat;

            % Plot spikes (purple)
            for spike = relative_spikes'
                plot([spike, spike], [y_pos - 0.4, y_pos + 0.4], '-', ...
                    'Color', color_us_catch, 'LineWidth', 0.5);
            end

            % US onset marker (dashed, outcome color) — placed at expected_us_lat
            % so it aligns vertically with normal trial outcome markers
            line([expected_us_lat, expected_us_lat], [y_pos - 0.4, y_pos + 0.4], ...
                'Color', outcome_color, 'LineStyle', '--', 'LineWidth', 1.5);
        end

        % =============================================================
        % Full-height vertical reference lines (behind spikes)
        % =============================================================

        % Odor 1 onset (red solid)
        line([0, 0], [0.5, max_y + 0.5], ...
            'Color', 'r', 'LineStyle', '-', 'LineWidth', 3, 'LineSmoothing', 'on');

        % Odor 2 onset at t = 1.5 s (blue dotted)
        line([1.5, 1.5], [0.5, max_y + 0.5], ...
            'Color', 'b', 'LineStyle', ':', 'LineWidth', 2.5, 'LineSmoothing', 'on');

        % Expected outcome time (outcome color, dashed)
        if ~isnan(expected_us_lat)
            line([expected_us_lat, expected_us_lat], [0.5, max_y + 0.5], ...
                'Color', outcome_color, 'LineStyle', '--', 'LineWidth', 2.5, 'LineSmoothing', 'on');
        end

        % Odor 1 offset at t = 1.0 s (gray)
        line([1.0, 1.0], [0.5, max_y + 0.5], ...
            'Color', [0.7 0.7 0.7], 'LineStyle', '-', 'LineWidth', 1.0);

        % Odor 2 offset at t = 2.5 s (gray)
        line([2.5, 2.5], [0.5, max_y + 0.5], ...
            'Color', [0.7 0.7 0.7], 'LineStyle', '-', 'LineWidth', 1.0);

        % =============================================================
        % Section separator lines and labels
        % =============================================================

        % Separator between normal and odor-only catch
        if n_normal > 0 && n_odor_catch > 0
            sep_y = (normal_y_end + odor_catch_y_start) / 2;
            line(xlim, [sep_y, sep_y], 'Color', [0.4 0.4 0.4], ...
                'LineStyle', '-', 'LineWidth', 1.5);
        end

        % Separator between odor-only catch and US-only catch
        if n_us_catch > 0
            if n_odor_catch > 0
                sep_y2 = (odor_catch_y_end + us_catch_y_start) / 2;
            else
                sep_y2 = (normal_y_end + us_catch_y_start) / 2;
            end
            line(xlim, [sep_y2, sep_y2], 'Color', [0.4 0.4 0.4], ...
                'LineStyle', '-', 'LineWidth', 1.5);
        end

        % Section labels on right side
        x_label_pos = max([reward_lat, airpuff_lat]) + time_after - 0.2;

        if n_normal > 0
            text(x_label_pos, (normal_y_start + normal_y_end) / 2, ...
                sprintf('Normal (n=%d)', n_normal), ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', color_normal, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end

        if n_odor_catch > 0
            text(x_label_pos, (odor_catch_y_start + odor_catch_y_end) / 2, ...
                sprintf('Odor-only catch (n=%d)', n_odor_catch), ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', color_odor_catch, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end

        if n_us_catch > 0
            text(x_label_pos, (us_catch_y_start + us_catch_y_end) / 2, ...
                sprintf('US-only catch (n=%d)', n_us_catch), ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', color_us_catch, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end

        % =============================================================
        % Legend
        % =============================================================
        legend_y = max_y + 1.5;
        legend_x_start = -time_before + 0.1;

        % Background box
        legend_width = 6.5;
        legend_height = 2.0;
        rectangle('Position', [legend_x_start - 0.3, legend_y - 1.0, legend_width, legend_height], ...
            'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'black', 'LineWidth', 1.5);

        % Row 1: Event markers
        % Odor 1 onset
        line([legend_x_start, legend_x_start + 0.2], [legend_y + 0.4, legend_y + 0.4], ...
            'Color', 'r', 'LineStyle', '-', 'LineWidth', 2.5);
        text(legend_x_start + 0.35, legend_y + 0.4, 'Odor 1 onset', 'FontSize', 9, ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold');

        % Odor 2 onset
        line([legend_x_start + 1.5, legend_x_start + 1.7], [legend_y + 0.4, legend_y + 0.4], ...
            'Color', 'b', 'LineStyle', ':', 'LineWidth', 2.5);
        text(legend_x_start + 1.85, legend_y + 0.4, 'Odor 2 onset', 'FontSize', 9, ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold');

        % Outcome
        line([legend_x_start + 3.0, legend_x_start + 3.2], [legend_y + 0.4, legend_y + 0.4], ...
            'Color', outcome_color, 'LineStyle', '--', 'LineWidth', 2.5);
        text(legend_x_start + 3.35, legend_y + 0.4, sprintf('%s', current_trial_type), ...
            'FontSize', 9, 'VerticalAlignment', 'middle', 'FontWeight', 'bold');

        % Row 2: Spike colors
        % Normal spikes
        line([legend_x_start, legend_x_start + 0.2], [legend_y - 0.2, legend_y - 0.2], ...
            'Color', color_normal, 'LineStyle', '-', 'LineWidth', 2.5);
        text(legend_x_start + 0.35, legend_y - 0.2, 'Normal trial', 'FontSize', 9, ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold');

        % Odor-only catch spikes
        line([legend_x_start + 1.5, legend_x_start + 1.7], [legend_y - 0.2, legend_y - 0.2], ...
            'Color', color_odor_catch, 'LineStyle', '-', 'LineWidth', 2.5);
        text(legend_x_start + 1.85, legend_y - 0.2, 'Odor-only catch', 'FontSize', 9, ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold');

        % US-only catch spikes
        line([legend_x_start + 3.5, legend_x_start + 3.7], [legend_y - 0.2, legend_y - 0.2], ...
            'Color', color_us_catch, 'LineStyle', '-', 'LineWidth', 2.5);
        text(legend_x_start + 3.85, legend_y - 0.2, 'US-only catch', 'FontSize', 9, ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold');

        % =============================================================
        % Axes formatting
        % =============================================================
        xlabel('Time relative to Odor 1 onset (s)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Trial Number', 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('Full Trial Raster - %s (%s trials)', ...
            replace(neuron_label, '_', ' '), current_trial_type), ...
            'FontSize', 14, 'FontWeight', 'bold');

        xlim([-time_before - 0.3, max([reward_lat, airpuff_lat]) + time_after + 0.3]);
        ylim([0, max_y + legend_height + 2.0]);

        grid on;
        grid minor;
        set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
        set(gca, 'FontSize', 11);

        hold off;

        % Save figure
        filename = fullfile(output_dir, sprintf('%s_full_trial_%s_raster_catch.png', ...
            neuron_label, current_trial_type));
        saveas(fig, filename);
        close(fig);

        fprintf('Saved: %s\n', filename);
    end
end

fprintf('All full trial raster plots generated and saved to %s\n', output_dir);
end
clear variables