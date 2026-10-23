%% SPIKE RASTER PLOT GENERATOR - FULL TRIAL VERSION
% This script generates raster plots for neural spiking data across full trials
% Each plot shows spike timing for a neuron across multiple trials of the same type
% (Reward trials, Airpuff trials, or Neutral trials)

%% Input information
cd('D:\vCA1\M656\preprocessed\M656_2026_03_25');
cfg.probe = 1;
cfg.spikeIDs = {'imec0_56'};
cfg.pre_trial = 1.0; % time before first odor to plot (seconds)
cfg.post_trial = 5.0; % time after outcome to plot (seconds) - INCREASED TO 5 SECONDS

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

%% Build trial structure using odor pair information
trials = struct('odor1_onset', {}, 'odor1_name', {}, 'odor2_onset', {}, 'odor2_name', {}, ...
                'outcome_time', {}, 'outcome_type', {});
trial_count = 0;

% Define odor pairs (based on your description: A→B, C→D, E→F)
odor_pairs = {};
if ismember('Odor A', [ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors]) && ismember('Odor B', [ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors])
    odor_pairs = [odor_pairs; {'Odor A', 'Odor B'}];
end
if ismember('Odor C', [ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors]) && ismember('Odor D', [ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors])
    odor_pairs = [odor_pairs; {'Odor C', 'Odor D'}];
end
if ismember('Odor E', [ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors]) && ismember('Odor F', [ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors])
    odor_pairs = [odor_pairs; {'Odor E', 'Odor F'}];
end

% Process each odor pair
for pair_idx = 1:size(odor_pairs, 1)
    odor1_name = odor_pairs{pair_idx, 1};
    odor2_name = odor_pairs{pair_idx, 2};
    
    % Find all onsets of odor1
    odor1_mask = strcmp(odor_identities, odor1_name);
    odor1_onsets = odor_onsets(odor1_mask);
    
    % For each odor1 onset, find the next odor2 onset
    for i = 1:length(odor1_onsets)
        odor1_time = odor1_onsets(i);
        
        % Find the next odor2 onset that occurs after this odor1
        subsequent_odor2 = odor_onsets(strcmp(odor_identities, odor2_name) & odor_onsets > odor1_time);
        
        if ~isempty(subsequent_odor2)
            % Use the first subsequent odor2
            odor2_time = subsequent_odor2(1);
            
            % Check that this is a valid pair (odor2 should be close in time)
            if odor2_time - odor1_time < 5
                trial_count = trial_count + 1;
                trials(trial_count).odor1_onset = odor1_time;
                trials(trial_count).odor1_name = odor1_name;
                trials(trial_count).odor2_onset = odor2_time;
                trials(trial_count).odor2_name = odor2_name;
                
                % Determine outcome type based on the second odor's category
                if isKey(odor_to_outcome, odor2_name)
                    outcome_type = odor_to_outcome(odor2_name);
                else
                    outcome_type = 'Unknown';
                end
                trials(trial_count).outcome_type = outcome_type;
                
                % Find the corresponding US time based on outcome type
                outcome_window_start = odor2_time;
                outcome_window_end = odor2_time + 3;
                
                % Assign outcome time based on predicted outcome type
                switch outcome_type
                    case 'Reward'
                        reward_in_window = reward_times(reward_times >= outcome_window_start & ...
                                                        reward_times <= outcome_window_end);
                        if ~isempty(reward_in_window)
                            trials(trial_count).outcome_time = reward_in_window(1);
                        else
                            trials(trial_count).outcome_time = outcome_window_start + 1;
                        end
                        
                    case 'Airpuff'
                        airpuff_in_window = airpuff_times(airpuff_times >= outcome_window_start & ...
                                                          airpuff_times <= outcome_window_end);
                        if ~isempty(airpuff_in_window)
                            trials(trial_count).outcome_time = airpuff_in_window(1);
                        else
                            trials(trial_count).outcome_time = outcome_window_start + 1;
                        end
                        
                    case 'Neutral'
                        trials(trial_count).outcome_time = outcome_window_start + 1;
                        
                    otherwise
                        trials(trial_count).outcome_time = outcome_window_start + 1;
                end
            end
        end
    end
end

fprintf('Found %d total trials\n', length(trials));
fprintf('  Reward trials: %d\n', sum(strcmp({trials.outcome_type}, 'Reward')));
fprintf('  Airpuff trials: %d\n', sum(strcmp({trials.outcome_type}, 'Airpuff')));
fprintf('  Neutral trials: %d\n', sum(strcmp({trials.outcome_type}, 'Neutral')));

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
    generate_full_trial_rasters(S_plot, trials, cfg.pre_trial, cfg.post_trial, './full_trial_rasters');
else
    warning('No valid trials found. Check event file and odor labels.');
end

%% Clear variables
% clearvars

%% Function to generate full trial raster plots
function generate_full_trial_rasters(S, trials, time_before, time_after, output_dir)
% GENERATE_FULL_TRIAL_RASTERS Creates raster plots for neurons across full trials

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Define trial types
trial_types = {'Reward', 'Airpuff', 'Neutral'};

% Define colors for different outcome types
outcome_colors = containers.Map({'Reward', 'Airpuff', 'Neutral'}, ...
                                {[0, 0.8, 0], [1, 0, 0], [0.5, 0.5, 0.5]});

% Loop through each neuron
for neuron_idx = 1:length(S.t)
    spike_times = S.t{neuron_idx};
    neuron_label = S.label{neuron_idx};
    
    % Create three figures: one for each trial type
    for trial_type_idx = 1:length(trial_types)
        current_trial_type = trial_types{trial_type_idx};
        
        % Filter trials by type
        trial_indices = find(strcmp({trials.outcome_type}, current_trial_type));
        
        if isempty(trial_indices)
            fprintf('No %s trials found for %s\n', current_trial_type, neuron_label);
            continue;
        end
        
        % Create figure with larger size for better visibility
        fig = figure('Position', [100, 100, 1600, 800]);
        hold on;
        
        % Get outcome color for this trial type
        outcome_color = outcome_colors(current_trial_type);
        
        % Plot each trial
        for plot_idx = 1:length(trial_indices)
            trial_num = trial_indices(plot_idx);
            trial = trials(trial_num);
            
            % Define trial window
            trial_start = trial.odor1_onset - time_before;
            trial_end = trial.outcome_time + time_after;
            
            % Find spikes within trial window
            spikes_in_trial = spike_times(spike_times >= trial_start & spike_times <= trial_end);
            
            % Convert to relative time (0 = odor1 onset)
            relative_spike_times = spikes_in_trial - trial.odor1_onset;
            
            % Plot spikes as vertical lines (black)
            for spike = relative_spike_times'
                plot([spike, spike], [plot_idx - 0.4, plot_idx + 0.4], 'k-', 'LineWidth', 0.5);
            end
            
            % Calculate times relative to odor1 onset
            odor2_relative = trial.odor2_onset - trial.odor1_onset;
            outcome_relative = trial.outcome_time - trial.odor1_onset;
            
            % Odor 1 onset marker (red solid line at x=0)
            line([0, 0], [plot_idx - 0.4, plot_idx + 0.4], ...
                'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
            
            % Odor 2 onset marker (blue dotted line)
            line([odor2_relative, odor2_relative], [plot_idx - 0.4, plot_idx + 0.4], ...
                'Color', 'b', 'LineStyle', ':', 'LineWidth', 1.5);
            
            % Outcome marker (dashed line, color depends on type)
            line([outcome_relative, outcome_relative], [plot_idx - 0.4, plot_idx + 0.4], ...
                'Color', outcome_color, 'LineStyle', '--', 'LineWidth', 1.5);
        end
        
        % Add full-height vertical reference lines behind spikes
        max_trial = length(trial_indices);
        
        % Get the reference trial to extract timing information
        ref_trial = trials(trial_indices(1));
        odor2_relative = ref_trial.odor2_onset - ref_trial.odor1_onset;
        outcome_relative = ref_trial.outcome_time - ref_trial.odor1_onset;
        
        % Odor 1 onset vertical line (red, solid, thicker)
        line([0, 0], [0.5, max_trial + 0.5], ...
            'Color', 'r', 'LineStyle', '-', 'LineWidth', 3, 'LineSmoothing', 'on');
        
        % Odor 2 onset vertical line (blue, dotted, thicker)
        line([odor2_relative, odor2_relative], [0.5, max_trial + 0.5], ...
            'Color', 'b', 'LineStyle', ':', 'LineWidth', 2.5, 'LineSmoothing', 'on');
        
        % Outcome vertical line (outcome color, dashed, thicker)
        line([outcome_relative, outcome_relative], [0.5, max_trial + 0.5], ...
            'Color', outcome_color, 'LineStyle', '--', 'LineWidth', 2.5, 'LineSmoothing', 'on');
        
        % Add legend box with key stimulus times
        legend_y = max_trial + 1.0;
        legend_x_start = -time_before + 0.1;
        
        % Create semi-transparent background box for legend
        rect = rectangle('Position', [legend_x_start - 0.3, legend_y - 0.6, 3.5, 1.2], ...
                         'FaceColor', [1 1 1], 'EdgeColor', 'black', 'LineWidth', 1.5);
        rect.FaceColor = [0.95 0.95 0.95];
        
        % Odor 1 onset legend
        line([legend_x_start, legend_x_start + 0.2], [legend_y, legend_y], ...
            'Color', 'r', 'LineStyle', '-', 'LineWidth', 2.5);
        text(legend_x_start + 0.35, legend_y, 'Odor 1 onset', 'FontSize', 11, ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold');
        
        % Odor 2 onset legend
        line([legend_x_start + 1.2, legend_x_start + 1.4], [legend_y - 0.4, legend_y - 0.4], ...
            'Color', 'b', 'LineStyle', ':', 'LineWidth', 2.5);
        text(legend_x_start + 1.55, legend_y - 0.4, 'Odor 2 onset', 'FontSize', 11, ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold');
        
        % Outcome legend
        line([legend_x_start + 2.4, legend_x_start + 2.6], [legend_y - 0.8, legend_y - 0.8], ...
            'Color', outcome_color, 'LineStyle', '--', 'LineWidth', 2.5);
        text(legend_x_start + 2.75, legend_y - 0.8, sprintf('%s', current_trial_type), ...
            'FontSize', 11, 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
        
        % Labels and formatting
        xlabel('Time relative to Odor 1 onset (s)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Trial Number', 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('Full Trial Raster - %s (%s trials, n=%d)', ...
            replace(neuron_label, '_', ' '), current_trial_type, length(trial_indices)), ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        xlim([-time_before - 0.3, time_after + 0.3]);
        ylim([0, max_trial + 1.8]);
        
        grid on;
        grid minor;
        set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
        set(gca, 'FontSize', 11);
        
        hold off;
        
        % Save figure as PNG with high resolution
        filename = fullfile(output_dir, sprintf('%s_full_trial_%s_raster.png', neuron_label, current_trial_type));
        saveas(fig, filename);
        close(fig);
        
        fprintf('Saved: %s\n', filename);
    end
end

fprintf('All full trial raster plots generated and saved to %s\n', output_dir);
end