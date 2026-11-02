%% Generates cell fingerprints for Odorpixels task 4 (odor, reward, airpuff) with catch trials
%  ADAPTED for catch trials + full-trial PETH
%
% OVERVIEW
% This script generates cell fingerprints for neurons recorded
% during Odorpixels task 4 with catch trials. Trial types include:
%   - Normal trials: odor -> outcome (reward/airpuff/dummyClick(neutral))
%   - Odor-only catch: odor -> NO outcome (i.e. dummyClick when reward or airpuff was expected)
%   - US-only catch: reward/airpuff delivered with NO preceding odor
%
% Each outcome category (reward, airpuff, neutral) has TWO associated odors
% (ExpKeys.reward_odors, ExpKeys.airpuff_odors, ExpKeys.neutral_odors).
% Either odor can start a trial, so onset times of the two odors are pooled
% within each category.
%
% dummyClick (evt label 'DummyClick') is delivered any time the reward or
% airpuff is NOT delivered: on neutral trials, and on odor-only catch trials.
%
% TRIAL TIMING (from odor onset):
%   - Odor onset:          0.00 s
%   - Odor offset:       2.00 s
%   - Airpuff/reward/dummyClick delivery: 1 s after odor offset (3.00 s from
%     odor onset), jittered by 200 ms (taken as +/- 0.2 s, i.e. 2.8-3.2 s)
%
% ODOR PETH (full trial)
%   Aligned to odor onset, window [-1, 4] s
%   One line PER ODOR (6 odors), colored by the outcome the odor predicts
%   (blue = reward odors, red = airpuff odors, black/gray = neutral odors).
%   Catch and paired trials of the same odor are pooled.
%
% US PETH
%   Aligned to outcome delivery, window [-2, 2] s
%   Reward, airpuff and dummyClick, each pooled over all deliveries
%   (paired + catch)
%
% DEPENDENCIES
% - Requires multiple functions from the mvdmlab and odor-pixels codebase
% - Data files:
%   * clean_units_imec0.mat (required)
%   * clean_units_imec1.mat (optional, for dual-probe recordings)
%   * imec0_clean_lfp.mat
%   * all_events.mat
%   * ExpKeys file (loaded via LoadExpKeys)
%
% OUTPUT
% PNG images of cell fingerprints saved with naming convention:
% [subject]_[date]_[unit_label]_cell_fingerprint.png

%% Input information (folder list)
folder_list = {'D:\Task4\M691\preprocessed\M691_2026_08_13'};

for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});
    % Loading data
    rng(2024);
    LoadExpKeys;

    % One probe obviously exists
    cfg_master = [];
    cfg_md = [];
    record_depth = ExpKeys.probe1_Depth;
    cfg_md.maxdepth = 1000 * record_depth;

    % Load spikes
    S = load_NPX_spike('clean_units_imec0.mat');
    S.usr.recording_depths = (cfg_md.maxdepth - S.usr.recording_depths) * cos(deg2rad(ExpKeys.probe1_insertion_roll));

    % If "imec1" exits, load the neurons from that as well
    if isfile('clean_units_imec1.mat')
        record_depth2 = ExpKeys.probe2_Depth;
        S2 = load_NPX_spike('clean_units_imec1.mat');
        S2.usr.recording_depths = (1000 * record_depth2 - S2.usr.recording_depths) * cos(deg2rad(ExpKeys.probe2_insertion_roll));
        S.label = cat(2, S.label, S2.label);
        S.t = cat(2, S.t, S2.t);
        S.usr.recording_depths = cat(2, S.usr.recording_depths, S2.usr.recording_depths);
        S.usr.shank_ids = cat(2, S.usr.shank_ids, S2.usr.shank_ids);
        S.usr.mean_waveforms = cat(2, S.usr.mean_waveforms, S2.usr.mean_waveforms);
        if (isfield(S.usr, 'channel_ids')) & (isfield(S2.usr, 'channel_ids'))
            S.usr.channel_ids = cat(2, S.usr.channel_ids, S2.usr.channel_ids);
        end
        clear S2
    end

    % Making binned version
    cfg_Q = [];
    cfg_Q.dt = 1;
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_size = 5;
    cfg_Q.gausswin_sd = 1;
    cfg_Q.boxcar_size = [];
    cfg_Q.tvec_edges = [];
    cfg_Q.verbose = 1;
    Q = MakeQfromS(cfg_Q, S);

    %%  Load LFPs and convert to mvdmlab format
    csc = tsd;
    load('imec0_clean_lfp.mat');

    csc.tvec = imec0.lfp_tvec';
    csc.data = double(imec0.lfp_traces');
    csc.label = cellstr(string(imec0.channel_ids).strip().extractAfter('ap#')');
    csc.cfg.depths = imec0.depths;

    Fs = struct;
    Fs.SamplingFrequency = double(imec0.lfp_fs);
    csc.cfg.hdr{1,1} = Fs;
    csc.cfg.depths = (cfg_md.maxdepth - csc.cfg.depths) * cos(deg2rad(ExpKeys.probe1_insertion_roll));

    % Order LFPs from most dorsal to ventral for sensible plotting
    [~, sort_idx] = sort(csc.cfg.depths);
    csc.data = csc.data(sort_idx,:);
    csc.label = csc.label(sort_idx);
    csc.cfg.depths = csc.cfg.depths(sort_idx);
    clear sort_idx

    % Calculate MUA for entire session
    cfg_MUA = [];
    cfg_MUA.tvec = csc.tvec;
    MUA = getMUA(cfg_MUA, S);
    MUAz = zscore_tsd(MUA);

    %%  Load events and classify what odors are what
    load('all_events.mat')

    % Create odor category mapping
    reward_indices = [];
    airpuff_indices = [];
    neutral_indices = [];

    if isfield(ExpKeys, 'reward_odors')
        for i = 1:length(evt.label)
            if any(strcmp(evt.label{i}, ExpKeys.reward_odors))
                reward_indices = [reward_indices, i];
            end
        end
    end

    if isfield(ExpKeys, 'airpuff_odors')
        for i = 1:length(evt.label)
            if any(strcmp(evt.label{i}, ExpKeys.airpuff_odors))
                airpuff_indices = [airpuff_indices, i];
            end
        end
    end

    if isfield(ExpKeys, 'neutral_odors')
        for i = 1:length(evt.label)
            if any(strcmp(evt.label{i}, ExpKeys.neutral_odors))
                neutral_indices = [neutral_indices, i];
            end
        end
    end

    %% Odor labels for each trial category (used in the figure title)
    % Each ExpKeys entry holds the TWO odors associated with that outcome
    % e.g. ExpKeys.reward_odors = {'Odor E', 'Odor F'} -> both are reward odors
    reward_odors_labels  = ExpKeys.reward_odors;
    reward_odors_labels = strjoin(extractAfter(reward_odors_labels, 'Odor '), '');
    airpuff_odors_labels = ExpKeys.airpuff_odors;
    airpuff_odors_labels = strjoin(extractAfter(airpuff_odors_labels, 'Odor '), '');
    neutral_odors_labels = ExpKeys.neutral_odors;
    neutral_odors_labels = strjoin(extractAfter(neutral_odors_labels, 'Odor '), '');

    %% Pool odor onset times for each trial category
    % reward_indices / airpuff_indices / neutral_indices (found above) are the
    % rows of evt belonging to the two odors of each category. Either odor starts
    % a trial, so concatenate their onset times, then sort so that trials are
    % in chronological order (the two odors are interleaved in the session).
    reward_odor_times = [];
    for i = reward_indices
        reward_odor_times = [reward_odor_times; evt.t{i}(:)];
    end
    reward_odor_times = sort(reward_odor_times);

    airpuff_odor_times = [];
    for i = airpuff_indices
        airpuff_odor_times = [airpuff_odor_times; evt.t{i}(:)];
    end
    airpuff_odor_times = sort(airpuff_odor_times);

    neutral_odor_times = [];
    for i = neutral_indices
        neutral_odor_times = [neutral_odor_times; evt.t{i}(:)];
    end
    neutral_odor_times = sort(neutral_odor_times);

    %% Load unconditioned stimulus events
    reward_times = evt.t(strcmp(evt.label, 'Reward'));
    airpuff_times = evt.t(strcmp(evt.label, 'Airpuff'));
    dummyClick_times = evt.t(strcmp(evt.label, 'DummyClick'));

    if ~isempty(reward_times)
        reward_times = reward_times{1};
    else
        reward_times = [];
    end

    if ~isempty(airpuff_times)
        airpuff_times = airpuff_times{1};
    else
        airpuff_times = [];
    end

    if ~isempty(dummyClick_times)
        dummyClick_times = dummyClick_times{1};
    else
        dummyClick_times = [];
    end

    %% Outcome timing (measured from odor onset)
    % 2 s odor + 1 s delay = 3 s, jittered by 200 ms (taken as +/- 0.2 s)
    outcome_latency_from_odor = 3.0;
    outcome_jitter = 0.2;

    %% Classify paired vs catch trials (diagnostic)
    % NOTE: this only affects the printed counts below. The PETHs pool all
    % trials of each category (paired + catch).
    % An outcome is "paired" with an odor if it falls within the jitter window
    % around outcome_latency_from_odor after that odor's onset. The window is
    % strictly positive, so outcomes BEFORE the odor can never match.
    pair_tol = outcome_jitter;

    % Classify reward USs: paired if a reward odor preceded at ~3.0 s
    reward_paired_times = [];
    reward_catch_times  = [];
    for i = 1:length(reward_times)
        dt = reward_times(i) - reward_odor_times;
        if any(abs(dt - outcome_latency_from_odor) <= pair_tol)
            reward_paired_times(end+1,1) = reward_times(i);
        else
            reward_catch_times(end+1,1)  = reward_times(i);
        end
    end

    % Classify airpuff USs
    airpuff_paired_times = [];
    airpuff_catch_times  = [];
    for i = 1:length(airpuff_times)
        dt = airpuff_times(i) - airpuff_odor_times;
        if any(abs(dt - outcome_latency_from_odor) <= pair_tol)
            airpuff_paired_times(end+1,1) = airpuff_times(i);
        else
            airpuff_catch_times(end+1,1)  = airpuff_times(i);
        end
    end

    % Classify dummyClicks by the odor that preceded them:
    %   neutral odor           -> expected neutral outcome
    %   reward/airpuff odor    -> odor-only catch (outcome omitted)
    %   no odor at ~3.0 s      -> unmatched
    dummyClick_neutral = 0; dummyClick_omission = 0; dummyClick_unmatched = 0;
    reward_airpuff_odor_times = [reward_odor_times; airpuff_odor_times];
    for i = 1:length(dummyClick_times)
        dt_neutral = dummyClick_times(i) - neutral_odor_times;
        dt_rewpuff = dummyClick_times(i) - reward_airpuff_odor_times;
        if any(abs(dt_neutral - outcome_latency_from_odor) <= pair_tol)
            dummyClick_neutral   = dummyClick_neutral + 1;
        elseif any(abs(dt_rewpuff - outcome_latency_from_odor) <= pair_tol)
            dummyClick_omission  = dummyClick_omission + 1;
        else
            dummyClick_unmatched = dummyClick_unmatched + 1;
        end
    end

    % Classify odor trials by what followed: expected outcome (paired),
    % dummyClick (odor-only catch), or neither (other, e.g. trial cut off at session end)
    reward_odor_paired = 0; reward_odor_catch = 0; reward_odor_other = 0;
    for i = 1:length(reward_odor_times)
        dt_us = reward_times - reward_odor_times(i);
        dt_dc = dummyClick_times - reward_odor_times(i);
        if any(abs(dt_us - outcome_latency_from_odor) <= pair_tol)
            reward_odor_paired = reward_odor_paired + 1;
        elseif any(abs(dt_dc - outcome_latency_from_odor) <= pair_tol)
            reward_odor_catch  = reward_odor_catch  + 1;
        else
            reward_odor_other  = reward_odor_other  + 1;
        end
    end
    airpuff_odor_paired = 0; airpuff_odor_catch = 0; airpuff_odor_other = 0;
    for i = 1:length(airpuff_odor_times)
        dt_us = airpuff_times - airpuff_odor_times(i);
        dt_dc = dummyClick_times - airpuff_odor_times(i);
        if any(abs(dt_us - outcome_latency_from_odor) <= pair_tol)
            airpuff_odor_paired = airpuff_odor_paired + 1;
        elseif any(abs(dt_dc - outcome_latency_from_odor) <= pair_tol)
            airpuff_odor_catch  = airpuff_odor_catch  + 1;
        else
            airpuff_odor_other  = airpuff_odor_other  + 1;
        end
    end
    neutral_odor_dummyClick = 0; neutral_odor_other = 0;
    for i = 1:length(neutral_odor_times)
        dt_dc = dummyClick_times - neutral_odor_times(i);
        if any(abs(dt_dc - outcome_latency_from_odor) <= pair_tol)
            neutral_odor_dummyClick = neutral_odor_dummyClick + 1;
        else
            neutral_odor_other      = neutral_odor_other + 1;
        end
    end

    fprintf('\n=== Catch trial diagnostics ===\n');
    fprintf('Reward  USs: %d total (%d paired, %d catch/US-only)\n', ...
        length(reward_times), length(reward_paired_times), length(reward_catch_times));
    fprintf('Airpuff USs: %d total (%d paired, %d catch/US-only)\n', ...
        length(airpuff_times), length(airpuff_paired_times), length(airpuff_catch_times));
    fprintf('DummyClicks: %d total (%d after neutral odor, %d odor-only catch, %d unmatched)\n', ...
        length(dummyClick_times), dummyClick_neutral, dummyClick_omission, dummyClick_unmatched);
    fprintf('Reward  odor trials:  %d total (%d paired, %d odor-only catch, %d other)\n', ...
        length(reward_odor_times), reward_odor_paired, reward_odor_catch, reward_odor_other);
    fprintf('Airpuff odor trials:  %d total (%d paired, %d odor-only catch, %d other)\n', ...
        length(airpuff_odor_times), airpuff_odor_paired, airpuff_odor_catch, airpuff_odor_other);
    fprintf('Neutral odor trials:  %d total (%d with dummyClick, %d other)\n', ...
        length(neutral_odor_times), neutral_odor_dummyClick, neutral_odor_other);
    fprintf('===============================\n\n');

    %%  Set up colors and PETH parameters

    % Full-trial PETH config (aligned to odor onset)
    cfg_peth_trial = [];
    cfg_peth_trial.window = [-1 4];
    cfg_peth_trial.dt = 0.01;
    cfg_peth_trial.mode = 'interp';

    % US PETH config (aligned to US onset)
    cfg_peth_us = [];
    cfg_peth_us.window = [-2 4];
    cfg_peth_us.dt = 0.01;
    cfg_peth_us.mode = 'interp';

    % Colors for trial categories
    col_reward  = [0 0 1];    % blue
    col_airpuff = [1 0 0];    % red
    col_neutral = [0 0 0];    % black

    %% Create odor event objects (aligned to odor onset, one per odor)
    % Rows of evt for the two odors of each outcome (found from ExpKeys above)
    odor_indices = [reward_indices, airpuff_indices, neutral_indices];
    num_odors = length(odor_indices);

    % Outcome category of each odor, used in the legend
    odor_category = [repmat({'Reward'},  1, length(reward_indices)), ...
                     repmat({'Airpuff'}, 1, length(airpuff_indices)), ...
                     repmat({'Neutral'}, 1, length(neutral_indices))];

    % One color per odor: the two odors of a category are two shades of that
    % category's color (dark = first odor in evt, light = second odor)
    cols_reward  = [col_reward;  0.4 0.6 1.0];
    cols_airpuff = [col_airpuff; 1.0 0.6 0.6];
    cols_neutral = [col_neutral; 0.5 0.5 0.5];
    odor_cols = [cols_reward(1:length(reward_indices), :); ...
                 cols_airpuff(1:length(airpuff_indices), :); ...
                 cols_neutral(1:length(neutral_indices), :)];

    odor_evts = cell(1, num_odors);
    for iE = 1:num_odors
        odor_evts{iE} = ts(evt.t(odor_indices(iE)), evt.label(odor_indices(iE)));

        % Remove events too close to recording edges
        if ~isempty(odor_evts{iE}.t{1})
            tstart = odor_evts{iE}.t{1} + cfg_peth_trial.window(1);
            tend   = odor_evts{iE}.t{1} + cfg_peth_trial.window(2);
            remove_idx = tstart < MUA.tvec(1) | tend > MUA.tvec(end);
            odor_evts{iE}.t{1}(remove_idx) = [];
        end
    end

    %% Create US event objects (all USs: paired + catch pooled)
    reward_evt  = ts({reward_times});
    airpuff_evt = ts({airpuff_times});
    dummyClick_evt = ts({dummyClick_times});

    % Remove US events too close to recording edges
    us_evts = {reward_evt, airpuff_evt, dummyClick_evt};
    for iE = 1:length(us_evts)
        if ~isempty(us_evts{iE}.t{1})
            tstart = us_evts{iE}.t{1} + cfg_peth_us.window(1);
            tend   = us_evts{iE}.t{1} + cfg_peth_us.window(2);
            remove_idx = tstart < MUA.tvec(1) | tend > MUA.tvec(end);
            us_evts{iE}.t{1}(remove_idx) = [];
        end
    end
    reward_evt  = us_evts{1};
    airpuff_evt = us_evts{2};
    dummyClick_evt = us_evts{3};

    %% Generate cell fingerprints
    for iC = 1:length(S.t)
        this_fig = figure('WindowState', 'maximized');
        this_cell = SelectTS([], S, iC);

        % Calculate firing rate and z-scored firing rate for this cell
        MUA_cell  = getMUA(cfg_MUA, this_cell);
        MUAz_cell = zscore_tsd(MUA_cell);

        % ============================================================
        % ODOR PETH (row 1, columns 1-4)
        % Aligned to odor onset, one line per odor
        % ============================================================
        ax1 = subplot(3, 5, [1, 2]);   % Raw full-trial PETH
        hold on
        ax2 = subplot(3, 5, [3, 4]);   % Z-scored full-trial PETH
        hold on

        leg_labels_trial = {};
        leg_cols_trial = [];
        min_ylim = 10;
        max_ylim = -10;

        % --- One line per odor ---
        for iOdor = 1:num_odors
            this_odor = restrict(odor_evts{iOdor}, iv(ExpKeys.sessstart, ExpKeys.sessend));
            if ~isempty(this_odor.t{1})
                out_z   = TSDpeth_fast(cfg_peth_trial, MUAz_cell, this_odor);
                out_raw = TSDpeth_fast(cfg_peth_trial, MUA_cell,  this_odor);
                leg_labels_trial{end+1} = sprintf('%s: %s (n=%d)', odor_category{iOdor}, ...
                    odor_evts{iOdor}.label{1}, length(this_odor.t{1}));
                leg_cols_trial(end+1, :) = odor_cols(iOdor, :);

                plot(ax1, out_raw.tvec, out_raw.data, 'Color', odor_cols(iOdor, :), 'LineWidth', 2);
                plot(ax2, out_z.tvec,   out_z.data,   'Color', odor_cols(iOdor, :), 'LineWidth', 2);
                min_ylim = min(min_ylim, min(out_z.data) - 0.1);
                max_ylim = max(max_ylim, max(out_z.data) + 0.1);
            end
        end

        % --- Reference lines for trial structure ---
        for ax_curr = [ax1, ax2]
            % Odor onset
            xline(ax_curr, 0,    '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', 'Odor on');
            % Odor offset
            xline(ax_curr, 2.0,  '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'Label', 'Odor off');
            % Expected outcome time (reward / airpuff / dummyClick all arrive at the same latency)
            xline(ax_curr, outcome_latency_from_odor, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', 'Outcome');
        end

        % Format raw full-trial PETH
        ax1.Title.String = 'Raw Odor PETH';
        ax1.Title.FontSize = 14;
        ax1.YLabel.String = 'F.R. (Hz)';
        ax1.YLabel.FontSize = 14;
        ax1.XLabel.String = 'Time from odor onset (s)';
        ax1.XLabel.FontSize = 14;

        % Format z-scored full-trial PETH
        ax2.YLim = [min(min_ylim, -0.5), max(max_ylim, 1)];
        ax2.Title.String = 'Z-scored Odor PETH';
        ax2.Title.FontSize = 14;
        ax2.YLabel.String = 'Z-score';
        ax2.YLabel.FontSize = 14;
        ax2.XLabel.String = 'Time from odor onset (s)';
        ax2.XLabel.FontSize = 14;

        % Legend for full-trial PETHs - place to the LEFT
        ax_legend_trial = axes('Position', [0.02, 0.65, 0.08, 0.25]);
        ax_legend_trial.Visible = 'off';

        h_dummy_trial = [];
        for iLeg = 1:length(leg_labels_trial)
            h_dummy_trial = [h_dummy_trial, line(NaN, NaN, 'Color', leg_cols_trial(iLeg, :), 'LineWidth', 2)];
        end
        legend(ax_legend_trial, h_dummy_trial, leg_labels_trial, 'Location', 'best', 'FontSize', 8, 'Box', 'off');

        % ============================================================
        % US PETH (row 2, columns 1-4) 
        % ============================================================
        ax3 = subplot(3, 5, [6, 7]);   % US Raw PETH
        hold on
        ax4 = subplot(3, 5, [8, 9]);   % US Z-scored PETH
        hold on

        leg_labels_us = {};
        min_ylim_us = 10;
        max_ylim_us = -10;

        % Plot reward US PETH
        reward_cue = restrict(reward_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(reward_cue.t{1})
            this_out_us_z   = TSDpeth_fast(cfg_peth_us, MUAz_cell, reward_cue);
            this_out_us_raw = TSDpeth_fast(cfg_peth_us, MUA_cell,  reward_cue);
            leg_labels_us{end+1} = sprintf('Reward (n=%d)', length(reward_cue.t{1}));

            plot(ax3, this_out_us_raw.tvec, this_out_us_raw.data, 'Color', col_reward, 'LineWidth', 2);
            plot(ax4, this_out_us_z.tvec,   this_out_us_z.data,   'Color', col_reward, 'LineWidth', 2);

            min_ylim_us = min(min_ylim_us, min(this_out_us_z.data) - 0.1);
            max_ylim_us = max(max_ylim_us, max(this_out_us_z.data) + 0.1);
        end

        % Plot airpuff US PETH
        airpuff_cue = restrict(airpuff_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(airpuff_cue.t{1})
            this_out_us_z   = TSDpeth_fast(cfg_peth_us, MUAz_cell, airpuff_cue);
            this_out_us_raw = TSDpeth_fast(cfg_peth_us, MUA_cell,  airpuff_cue);
            leg_labels_us{end+1} = sprintf('Airpuff (n=%d)', length(airpuff_cue.t{1}));

            plot(ax3, this_out_us_raw.tvec, this_out_us_raw.data, 'Color', col_airpuff, 'LineWidth', 2);
            plot(ax4, this_out_us_z.tvec,   this_out_us_z.data,   'Color', col_airpuff, 'LineWidth', 2);

            min_ylim_us = min(min_ylim_us, min(this_out_us_z.data) - 0.1);
            max_ylim_us = max(max_ylim_us, max(this_out_us_z.data) + 0.1);
        end

        % Plot dummyClick US PETH
        dummyClick_cue = restrict(dummyClick_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(dummyClick_cue.t{1})
            this_out_us_z   = TSDpeth_fast(cfg_peth_us, MUAz_cell, dummyClick_cue);
            this_out_us_raw = TSDpeth_fast(cfg_peth_us, MUA_cell,  dummyClick_cue);
            leg_labels_us{end+1} = sprintf('DummyClick (n=%d)', length(dummyClick_cue.t{1}));

            plot(ax3, this_out_us_raw.tvec, this_out_us_raw.data, 'Color', col_neutral, 'LineWidth', 2);
            plot(ax4, this_out_us_z.tvec,   this_out_us_z.data,   'Color', col_neutral, 'LineWidth', 2);

            min_ylim_us = min(min_ylim_us, min(this_out_us_z.data) - 0.1);
            max_ylim_us = max(max_ylim_us, max(this_out_us_z.data) + 0.1);
        end

        % Format US raw PETH plot
        xline(ax3, 0, ':', 'Color', [0 0 0], 'LineWidth', 1.5);
        ax3.Title.String = 'Raw US PETH';
        ax3.Title.FontSize = 14;
        ax3.YLabel.String = 'F.R. (Hz)';
        ax3.YLabel.FontSize = 14;
        ax3.XLabel.String = 'Time from US onset (s)';
        ax3.XLabel.FontSize = 14;

        % Format US z-scored PETH plot
        ax4.YLim = [min(min_ylim_us, -0.5), max(max_ylim_us, 1)];
        xline(ax4, 0, ':', 'Color', [0 0 0], 'LineWidth', 1.5);
        ax4.Title.String = 'Z-scored US PETH';
        ax4.Title.FontSize = 14;
        ax4.YLabel.String = 'Z-score';
        ax4.YLabel.FontSize = 14;
        ax4.XLabel.String = 'Time from US onset (s)';
        ax4.XLabel.FontSize = 14;

        % Legend for US PETHs - place to the LEFT
        ax_legend_us = axes('Position', [0.02, 0.35, 0.08, 0.25]);
        ax_legend_us.Visible = 'off';

        h_dummy_us = [];
        for iUS = 1:length(leg_labels_us)
            if contains(leg_labels_us{iUS}, 'Reward')
                h_dummy_us = [h_dummy_us, line(NaN, NaN, 'Color', col_reward, 'LineWidth', 2)];
            elseif contains(leg_labels_us{iUS}, 'Airpuff')
                h_dummy_us = [h_dummy_us, line(NaN, NaN, 'Color', col_airpuff, 'LineWidth', 2)];
            elseif contains(leg_labels_us{iUS}, 'DummyClick')
                h_dummy_us = [h_dummy_us, line(NaN, NaN, 'Color', col_neutral, 'LineWidth', 2)];
            end
        end
        legend(ax_legend_us, h_dummy_us, leg_labels_us, 'Location', 'best', 'FontSize', 8, 'Box', 'off');

        % ============================================================
        % WAVEFORM, AUTOCORRELATION, ISI, FIRING RATE
        % ============================================================

        % Plot waveform (row 1, column 5)
        subplot(3, 5, 5)
        [~, max_idx] = max((max(this_cell.usr.mean_waveforms{1}, [], 2)));
        plot(this_cell.usr.mean_waveforms{1}(max_idx,:))
        title('Mean Waveform', 'FontSize', 14)
        xticks([])
        yticks([])

        % Plot autocorrelation (row 2, column 5)
        subplot(3, 5, 10)
        cfg_xc = [];
        cfg_xc.smooth = 1;
        cfg_xc.max_t = 1;
        [test_xc, test_t] = ccf(cfg_xc, this_cell.t{1}, this_cell.t{1});
        plot(test_t, test_xc);
        xticks(0:1/8:1);
        xlim([0 1])
        xlabel('Time (s)', 'FontSize', 14)
        title('Auto-correlation', 'FontSize', 14)

        % Plot log ISI (row 3, column 5)
        subplot(3, 5, 15)
        test_isi = diff(this_cell.t{1});
        histogram(log10(test_isi), -3.75:0.5:2.75);
        xline(log10(0.002), 'r');
        this_str = sprintf(', %d/%d ISIs less than 2 msec', ...
            sum(test_isi < 0.002), length(test_isi));
        title(strcat('LogISI', this_str), 'FontSize', 14)
        ax = gca;
        ax.YAxis.Exponent = 0;
        xlabel('Time (ms)', 'FontSize', 14)

        % Plot firing-rate over time (row 3, columns 1-4)
        subplot(3, 5, [11:14])
        hold on
        plot(Q.tvec, Q.data(iC,:), 'black');
        xlabel('Time (s)', 'FontSize', 14)
        ylabel('F.R. (Hz)', 'FontSize', 14)
        xlim([0, ExpKeys.sessend])
        title('Firing Rate Over Time', 'FontSize', 14)

        fname = strcat(ExpKeys.subject, '_', strrep(ExpKeys.date, '-', '_'), ...
            '_shank', string(this_cell.usr.shank_ids), '_', this_cell.label{1}, '_cell_fingerprint_odorUS');
        tname = strcat(fname, 'Reward:', reward_odors_labels, " Airpuff:", airpuff_odors_labels, " Neutral:", neutral_odors_labels);
        sgtitle(tname, 'Interpreter', 'None');
        print(this_fig, '-dpng', fname)
        close;
    end

    clearvars -except folder_list iFolder
end