%% Classify trials for Odorpixels task 4 (odor, reward, airpuff) with catch trials
%
% Loads one session and labels every trial by:
%   - odor presented ('none' if no odor)
%   - outcome the odor is associated with (reward / airpuff / neutral)
%   - outcome that actually happened (reward / airpuff / neutral / none)
%   - trial type (normal / odor_catch / outcome_catch / unclassified)
%   - trial start and end time
%
% Trial timing follows cell_fingerprint_task4_odors_catch.m:
%   odor on (0 s) -> odor off (2 s) -> outcome at 3 s (+/- 0.2 s jitter)
%
% Trial types:
%   normal         odor followed by its own outcome (neutral odor -> dummy click)
%   odor_catch     reward/airpuff odor followed by a dummy click (outcome omitted)
%   outcome_catch  reward/airpuff delivered with no odor before it
%   unclassified   anything that fits none of the above
%
% The 'neutral' outcome is the dummy click.
%
% OUTPUT
% 'trials' is a struct array with one element per trial, in chronological order.
% It is saved to [subject]_[date]_trials.mat in the session folder.
%
% Example use in an analysis:
%   load('M691_2026_08_13_trials.mat');
%   is_normal = strcmp({trials.trial_type}, 'normal');
%   normal_trials = trials(is_normal);
%   normal_starts = [normal_trials.trial_start];


%% Input
session_folder = 'D:\Task4\M671\preprocessed\M671_2026_08_27';

%% Trial timing (seconds from odor onset)
outcome_latency = 3.0;      % 2 s odor + 1 s delay
outcome_jitter = 0.25;       % outcome can arrive within +/- 0.2 s of the latency w a 0.05s buffer
post_outcome_time = 0;      % extra time to include after the outcome in trial_end


%% Load session info and events
cd(session_folder);
LoadExpKeys;
load('all_events.mat');


%% Collect all odor onsets
% Each outcome has two odors. ExpKeys says which odors go with which outcome.
outcome_names = {'reward', 'airpuff', 'neutral'};
odors_for_outcome = {ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors};

odor_times = [];
odor_labels = {};
odor_outcomes = {};

for iOutcome = 1:length(outcome_names)
    this_odor_list = cellstr(odors_for_outcome{iOutcome});

    for iOdor = 1:length(this_odor_list)
        row = find(strcmp(evt.label, this_odor_list{iOdor}), 1);

        if ~isempty(row)
            these_times = evt.t{row}(:);

            % Keep only events inside the session
            these_times = these_times(these_times >= ExpKeys.sessstart & these_times <= ExpKeys.sessend);

            for iTime = 1:length(these_times)
                odor_times(end+1, 1) = these_times(iTime);
                odor_labels{end+1, 1} = this_odor_list{iOdor};
                odor_outcomes{end+1, 1} = outcome_names{iOutcome};
            end
        end
    end
end

% Put odors in chronological order (the odors are interleaved in the session)
[odor_times, sort_idx] = sort(odor_times);
odor_labels = odor_labels(sort_idx);
odor_outcomes = odor_outcomes(sort_idx);


%% Collect all outcome times (reward, airpuff, dummy click)
% These are the event labels in evt, in the same order as outcome_names
outcome_labels_in_evt = {'Reward', 'Airpuff', 'DummyClick'};

outcome_times = [];
outcome_types = {};

for iOutcome = 1:length(outcome_labels_in_evt)
    row = find(strcmp(evt.label, outcome_labels_in_evt{iOutcome}), 1);

    if ~isempty(row)
        these_times = evt.t{row}(:);
        these_times = these_times(these_times >= ExpKeys.sessstart & these_times <= ExpKeys.sessend);

        for iTime = 1:length(these_times)
            outcome_times(end+1, 1) = these_times(iTime);
            outcome_types{end+1, 1} = outcome_names{iOutcome};
        end
    end
end

[outcome_times, sort_idx] = sort(outcome_times);
outcome_types = outcome_types(sort_idx);


%% Classify trials that start with an odor
trials = struct('trial_num', {}, 'odor', {}, 'odor_outcome', {}, 'actual_outcome', {}, ...
    'trial_type', {}, 'trial_start', {}, 'trial_end', {}, 'outcome_time', {});
n_trials = 0;

% Keeps track of which outcomes were paired with an odor
outcome_used = false(length(outcome_times), 1);

for iTrial = 1:length(odor_times)
    odor_on = odor_times(iTrial);
    expected_outcome = odor_outcomes{iTrial};

    % Look for an outcome about 3 s after this odor onset
    timing_diff = abs(outcome_times - odor_on - outcome_latency);
    match_idx = find(timing_diff <= outcome_jitter, 1);

    if isempty(match_idx)
        % No outcome arrived. Trial end uses the time the outcome was expected
        actual_outcome = 'none';
        outcome_time = NaN;
        trial_end = odor_on + outcome_latency + post_outcome_time;
    else
        actual_outcome = outcome_types{match_idx};
        outcome_time = outcome_times(match_idx);
        trial_end = outcome_time + post_outcome_time;
        outcome_used(match_idx) = true;
    end

    % Decide trial type
    if strcmp(actual_outcome, expected_outcome)
        trial_type = 'normal';
    elseif strcmp(actual_outcome, 'neutral') && ~strcmp(expected_outcome, 'neutral')
        trial_type = 'odor_catch';
    else
        trial_type = 'unclassified';
    end

    n_trials = n_trials + 1;
    trials(n_trials).odor = odor_labels{iTrial};
    trials(n_trials).odor_outcome = expected_outcome;
    trials(n_trials).actual_outcome = actual_outcome;
    trials(n_trials).trial_type = trial_type;
    trials(n_trials).trial_start = odor_on;
    trials(n_trials).trial_end = trial_end;
    trials(n_trials).outcome_time = outcome_time;
end


%% Classify outcomes with no odor before them
% Reward/airpuff with no odor = outcome catch. There is no odor onset, so the
% trial start is set to where the odor would have started in a normal trial.
for iOutcome = 1:length(outcome_times)
    if ~outcome_used(iOutcome)

        if strcmp(outcome_types{iOutcome}, 'neutral')
            % A dummy click with no odor is not an expected trial
            trial_type = 'unclassified';
        else
            trial_type = 'outcome_catch';
        end

        outcome_time = outcome_times(iOutcome);

        n_trials = n_trials + 1;
        trials(n_trials).odor = 'none';
        trials(n_trials).odor_outcome = 'none';
        trials(n_trials).actual_outcome = outcome_types{iOutcome};
        trials(n_trials).trial_type = trial_type;
        trials(n_trials).trial_start = outcome_time - outcome_latency;
        trials(n_trials).trial_end = outcome_time + post_outcome_time;
        trials(n_trials).outcome_time = outcome_time;
    end
end


%% Sort all trials by start time and number them
[~, sort_idx] = sort([trials.trial_start]);
trials = trials(sort_idx);

for iTrial = 1:length(trials)
    trials(iTrial).trial_num = iTrial;
end


%% Save
save_name = 'all_trials.mat';
save(fullfile(session_folder, save_name), 'trials');


%% Print trial counts
trial_type_list = {'normal', 'odor_catch', 'outcome_catch', 'unclassified'};
all_trial_types = {trials.trial_type};

fprintf('\n=== Trial counts: %s ===\n', save_name);
for iType = 1:length(trial_type_list)
    n_this_type = sum(strcmp(all_trial_types, trial_type_list{iType}));
    fprintf('%-14s %d\n', trial_type_list{iType}, n_this_type);
end
fprintf('%-14s %d\n', 'total', length(trials));
