%% Generates cell fingerprints for Odorpixels task 3 (odor, reward, airpuff) with catch trials
%  ADAPTED for catch trials + full-trial PETH
%
% OVERVIEW
% This script generates cell fingerprints for neurons recorded
% during Odorpixels task 3 with catch trials. Trial types include:
%   - Normal trials: odor1 -> odor2 -> outcome (reward/airpuff/none)
%   - Odor-only catch: odor1 -> odor2 -> NO outcome (when one was expected)
%   - US-only catch: reward/airpuff delivered with NO preceding odor pair
%
% TRIAL TIMING (from odor 1 onset):
%   - Odor 1 onset:        0.00 s
%   - Odor 1 offset:       1.00 s
%   - Break:               1.00 - 1.50 s
%   - Odor 2 onset:        1.50 s
%   - Odor 2 offset:       2.50 s
%   - Airpuff delivery:    2.75 s (0.25 s after odor 2 offset)
%   - Reward delivery:     3.00 s (0.50 s after odor 2 offset)
%
% FULL-TRIAL PETH
%   Aligned to odor 1 onset, window [-1, 4] s
%   Three lines grouped by outcome category:
%     - Reward trials (blue)
%     - Airpuff trials (red)
%     - Neutral trials (black)
%   Odor 1 is identified as the FIRST odor listed in each ExpKeys pair:
%     e.g. ExpKeys.reward_odors = {'Odor E', 'Odor F'} -> Odor E is odor 1
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
folder_list = {'D:\vCA1\M656\preprocessed\M656_2026_03_26_conservativeSpikeSorting'};

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

    %% Identify odor 1 labels for each trial category
    % The FIRST odor in each ExpKeys pair is odor 1 (the trial-initiating odor)
    % e.g. ExpKeys.reward_odors = {'Odor E', 'Odor F'} -> 'Odor E' is odor 1
    reward_odor1_label  = ExpKeys.reward_odors{1};
    airpuff_odor1_label = ExpKeys.airpuff_odors{1};
    neutral_odor1_label = ExpKeys.neutral_odors{1};

    reward_odors_labels  = ExpKeys.reward_odors;
    reward_odors_labels = strjoin(extractAfter(reward_odors_labels, 'Odor '), '');
    airpuff_odors_labels = ExpKeys.airpuff_odors;
    airpuff_odors_labels = strjoin(extractAfter(airpuff_odors_labels, 'Odor '), '');
    neutral_odors_labels = ExpKeys.neutral_odors;
    neutral_odors_labels = strjoin(extractAfter(neutral_odors_labels, 'Odor '), '');

    % Find the index in evt for each odor-1 label
    reward_odor1_idx  = find(strcmp(evt.label, reward_odor1_label));
    airpuff_odor1_idx = find(strcmp(evt.label, airpuff_odor1_label));
    neutral_odor1_idx = find(strcmp(evt.label, neutral_odor1_label));

    % Extract odor 1 onset times for each category
    if ~isempty(reward_odor1_idx)
        reward_odor1_times = evt.t{reward_odor1_idx};
    else
        reward_odor1_times = [];
    end

    if ~isempty(airpuff_odor1_idx)
        airpuff_odor1_times = evt.t{airpuff_odor1_idx};
    else
        airpuff_odor1_times = [];
    end

    if ~isempty(neutral_odor1_idx)
        neutral_odor1_times = evt.t{neutral_odor1_idx};
    else
        neutral_odor1_times = [];
    end

    %% Load unconditioned stimulus events
    reward_times = evt.t(strcmp(evt.label, 'Reward ON'));
    airpuff_times = evt.t(strcmp(evt.label, 'Airpuff ON'));

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

    %% Fixed US latencies (measured from odor 1 onset)
    % Odor 2 offset = 2.5 s from odor 1 onset
    % Airpuff = 2.5 + 0.25 = 2.75 s from odor 1 onset
    % Reward  = 2.5 + 0.50 = 3.0 s from odor 1 onset
    airpuff_latency_from_odor1 = 2.75;
    reward_latency_from_odor1  = 3.0;

    %% Classify paired vs catch trials (diagnostic)
    pair_tol = 0.5;  % tolerance in seconds around expected US time

    % Classify reward USs: paired if a reward odor-1 preceded at ~3.0 s
    reward_paired_times = [];
    reward_catch_times  = [];
    for i = 1:length(reward_times)
        dt = reward_times(i) - reward_odor1_times;
        dt = dt(dt > 0);
        if ~isempty(dt) && min(abs(dt - reward_latency_from_odor1)) <= pair_tol
            reward_paired_times(end+1,1) = reward_times(i);
        else
            reward_catch_times(end+1,1)  = reward_times(i);
        end
    end

    % Classify airpuff USs
    airpuff_paired_times = [];
    airpuff_catch_times  = [];
    for i = 1:length(airpuff_times)
        dt = airpuff_times(i) - airpuff_odor1_times;
        dt = dt(dt > 0);
        if ~isempty(dt) && min(abs(dt - airpuff_latency_from_odor1)) <= pair_tol
            airpuff_paired_times(end+1,1) = airpuff_times(i);
        else
            airpuff_catch_times(end+1,1)  = airpuff_times(i);
        end
    end

    % Classify odor trials: paired vs catch (odor-only)
    reward_odor_paired = 0; reward_odor_catch = 0;
    for i = 1:length(reward_odor1_times)
        dt = reward_times - reward_odor1_times(i);
        dt = dt(dt > 0);
        if ~isempty(dt) && min(abs(dt - reward_latency_from_odor1)) <= pair_tol
            reward_odor_paired = reward_odor_paired + 1;
        else
            reward_odor_catch  = reward_odor_catch  + 1;
        end
    end
    airpuff_odor_paired = 0; airpuff_odor_catch = 0;
    for i = 1:length(airpuff_odor1_times)
        dt = airpuff_times - airpuff_odor1_times(i);
        dt = dt(dt > 0);
        if ~isempty(dt) && min(abs(dt - airpuff_latency_from_odor1)) <= pair_tol
            airpuff_odor_paired = airpuff_odor_paired + 1;
        else
            airpuff_odor_catch  = airpuff_odor_catch  + 1;
        end
    end

    fprintf('\n=== Catch trial diagnostics ===\n');
    fprintf('Reward  USs: %d total (%d paired, %d catch/US-only)\n', ...
        length(reward_times), length(reward_paired_times), length(reward_catch_times));
    fprintf('Airpuff USs: %d total (%d paired, %d catch/US-only)\n', ...
        length(airpuff_times), length(airpuff_paired_times), length(airpuff_catch_times));
    fprintf('Reward  odor trials:  %d paired, %d odor-only catch\n', ...
        reward_odor_paired, reward_odor_catch);
    fprintf('Airpuff odor trials:  %d paired, %d odor-only catch\n', ...
        airpuff_odor_paired, airpuff_odor_catch);
    fprintf('Neutral odor trials:  %d total\n', length(neutral_odor1_times));
    fprintf('===============================\n\n');

    %%  Set up colors and PETH parameters

    % Full-trial PETH config (aligned to odor 1 onset)
    cfg_peth_trial = [];
    cfg_peth_trial.window = [-1 4];
    cfg_peth_trial.dt = 0.01;
    cfg_peth_trial.mode = 'interp';

    % US PETH config (aligned to US onset, same as before)
    cfg_peth_us = [];
    cfg_peth_us.window = [-2 2];
    cfg_peth_us.dt = 0.01;
    cfg_peth_us.mode = 'interp';

    % Colors for trial categories
    col_reward  = [0 0 1];    % blue
    col_airpuff = [1 0 0];    % red
    col_neutral = [0 0 0];    % black

    %% Create full-trial event objects (aligned to odor 1 onset)
    reward_trial_evt  = ts({reward_odor1_times});
    airpuff_trial_evt = ts({airpuff_odor1_times});
    neutral_trial_evt = ts({neutral_odor1_times});

    % Remove events too close to recording edges
    trial_evts = {reward_trial_evt, airpuff_trial_evt, neutral_trial_evt};
    for iE = 1:length(trial_evts)
        if ~isempty(trial_evts{iE}.t{1})
            tstart = trial_evts{iE}.t{1} + cfg_peth_trial.window(1);
            tend   = trial_evts{iE}.t{1} + cfg_peth_trial.window(2);
            remove_idx = tstart < MUA.tvec(1) | tend > MUA.tvec(end);
            trial_evts{iE}.t{1}(remove_idx) = [];
        end
    end
    reward_trial_evt  = trial_evts{1};
    airpuff_trial_evt = trial_evts{2};
    neutral_trial_evt = trial_evts{3};

    %% Create US event objects (all USs: paired + catch pooled)
    reward_evt  = ts({reward_times});
    airpuff_evt = ts({airpuff_times});

    % Remove US events too close to recording edges
    us_evts = {reward_evt, airpuff_evt};
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

    %% Generate cell fingerprints
    for iC = 1:length(S.t)
        this_fig = figure('WindowState', 'maximized');
        this_cell = SelectTS([], S, iC);

        % Calculate firing rate and z-scored firing rate for this cell
        MUA_cell  = getMUA(cfg_MUA, this_cell);
        MUAz_cell = zscore_tsd(MUA_cell);

        % ============================================================
        % FULL-TRIAL PETH (row 1, columns 1-4)
        % Aligned to odor 1 onset, grouped by outcome category
        % ============================================================
        ax1 = subplot(3, 5, [1, 2]);   % Raw full-trial PETH
        hold on
        ax2 = subplot(3, 5, [3, 4]);   % Z-scored full-trial PETH
        hold on

        leg_labels_trial = {};
        min_ylim = 10;
        max_ylim = -10;

        % --- Reward trials ---
        reward_trial_restricted = restrict(reward_trial_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(reward_trial_restricted.t{1})
            out_z   = TSDpeth_fast(cfg_peth_trial, MUAz_cell, reward_trial_restricted);
            out_raw = TSDpeth_fast(cfg_peth_trial, MUA_cell,  reward_trial_restricted);
            leg_labels_trial{end+1} = sprintf('Reward (n=%d)', length(reward_trial_restricted.t{1}));

            plot(ax1, out_raw.tvec, out_raw.data, 'Color', col_reward, 'LineWidth', 2);
            plot(ax2, out_z.tvec,   out_z.data,   'Color', col_reward, 'LineWidth', 2);
            min_ylim = min(min_ylim, min(out_z.data) - 0.1);
            max_ylim = max(max_ylim, max(out_z.data) + 0.1);
        end

        % --- Airpuff trials ---
        airpuff_trial_restricted = restrict(airpuff_trial_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(airpuff_trial_restricted.t{1})
            out_z   = TSDpeth_fast(cfg_peth_trial, MUAz_cell, airpuff_trial_restricted);
            out_raw = TSDpeth_fast(cfg_peth_trial, MUA_cell,  airpuff_trial_restricted);
            leg_labels_trial{end+1} = sprintf('Airpuff (n=%d)', length(airpuff_trial_restricted.t{1}));

            plot(ax1, out_raw.tvec, out_raw.data, 'Color', col_airpuff, 'LineWidth', 2);
            plot(ax2, out_z.tvec,   out_z.data,   'Color', col_airpuff, 'LineWidth', 2);
            min_ylim = min(min_ylim, min(out_z.data) - 0.1);
            max_ylim = max(max_ylim, max(out_z.data) + 0.1);
        end

        % --- Neutral trials ---
        neutral_trial_restricted = restrict(neutral_trial_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(neutral_trial_restricted.t{1})
            out_z   = TSDpeth_fast(cfg_peth_trial, MUAz_cell, neutral_trial_restricted);
            out_raw = TSDpeth_fast(cfg_peth_trial, MUA_cell,  neutral_trial_restricted);
            leg_labels_trial{end+1} = sprintf('Neutral (n=%d)', length(neutral_trial_restricted.t{1}));

            plot(ax1, out_raw.tvec, out_raw.data, 'Color', col_neutral, 'LineWidth', 2);
            plot(ax2, out_z.tvec,   out_z.data,   'Color', col_neutral, 'LineWidth', 2);
            min_ylim = min(min_ylim, min(out_z.data) - 0.1);
            max_ylim = max(max_ylim, max(out_z.data) + 0.1);
        end

        % --- Reference lines for trial structure ---
        for ax_curr = [ax1, ax2]
            % Odor 1 onset
            xline(ax_curr, 0,    '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', 'Odor 1 on');
            % Odor 1 offset
            xline(ax_curr, 1.0,  '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'Label', 'Odor 1 off');
            % Odor 2 onset
            xline(ax_curr, 1.5,  '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Label', 'Odor 2 on');
            % Odor 2 offset
            xline(ax_curr, 2.5,  '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'Label', 'Odor 2 off');
            % Expected airpuff time
            xline(ax_curr, airpuff_latency_from_odor1, '--', 'Color', col_airpuff, 'LineWidth', 1.5, 'Label', 'Airpuff');
            % Expected reward time
            xline(ax_curr, reward_latency_from_odor1,  '--', 'Color', col_reward,  'LineWidth', 1.5, 'Label', 'Reward');
        end

        % Format raw full-trial PETH
        ax1.Title.String = 'Raw Full-Trial PETH';
        ax1.Title.FontSize = 14;
        ax1.YLabel.String = 'F.R. (Hz)';
        ax1.YLabel.FontSize = 14;
        ax1.XLabel.String = 'Time from Odor 1 onset (s)';
        ax1.XLabel.FontSize = 14;

        % Format z-scored full-trial PETH
        ax2.YLim = [min(min_ylim, -0.5), max(max_ylim, 1)];
        ax2.Title.String = 'Z-scored Full-Trial PETH';
        ax2.Title.FontSize = 14;
        ax2.YLabel.String = 'Z-score';
        ax2.YLabel.FontSize = 14;
        ax2.XLabel.String = 'Time from Odor 1 onset (s)';
        ax2.XLabel.FontSize = 14;

        % Legend for full-trial PETHs - place to the LEFT
        ax_legend_trial = axes('Position', [0.02, 0.65, 0.08, 0.25]);
        ax_legend_trial.Visible = 'off';

        h_dummy_trial = [];
        trial_cols = {col_reward, col_airpuff, col_neutral};
        for iLeg = 1:length(leg_labels_trial)
            if contains(leg_labels_trial{iLeg}, 'Reward')
                h_dummy_trial = [h_dummy_trial, line(NaN, NaN, 'Color', col_reward, 'LineWidth', 2)];
            elseif contains(leg_labels_trial{iLeg}, 'Airpuff')
                h_dummy_trial = [h_dummy_trial, line(NaN, NaN, 'Color', col_airpuff, 'LineWidth', 2)];
            elseif contains(leg_labels_trial{iLeg}, 'Neutral')
                h_dummy_trial = [h_dummy_trial, line(NaN, NaN, 'Color', col_neutral, 'LineWidth', 2)];
            end
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
            '_', this_cell.label{1}, '_cell_fingerprint_odorUS');
        tname = strcat(fname, 'Reward:', reward_odors_labels, " Airpuff:", airpuff_odors_labels, " Neutral:", neutral_odors_labels);
        sgtitle(tname, 'Interpreter', 'None');
        print(this_fig, '-dpng', fname)
        close;
    end

    clearvars -except folder_list iFolder
end