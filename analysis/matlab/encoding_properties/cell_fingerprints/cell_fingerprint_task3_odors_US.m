%% Generates cell fingerprints for Odorpixels task 3 (odor, reward, airpuff)
%
% OVERVIEW
% This script generates cell fingerprints for neurons recorded
% during Odorpixels task 3 (no blocks). The cell fingerprint includes a
% visualization of single-unit neural activity, including PETHs, waveform
% characteristics, and firing rate dynamics.
%
% FUNCTIONALITY
% 1. Data Loading and Preprocessing: Loads spike data from cleaned Neuropixels
%    recordings (supporting single or dual-probe configurations), extracts local
%    field potentials (LFPs), and converts them to mvdmlab format for analysis.
%
% 2. Event Classification: Parses behavioral events from an events file and
%    classifies odor stimuli into three categories based on ExpKeys definitions:
%    - Reward-paired odors
%    - Airpuff-paired odors
%    - Neutral odors
%
% 3. Unconditioned Stimulus (US) Latency Calculation: Computes the mean latency
%    between the second odor in a pair and subsequent unconditioned stimuli
%    (reward or airpuff).
%
% 4. Cell Fingerprint Generation: For each neuron, creates a comprehensive
%    figure displaying:
%    - Z-scored and raw odor PETHs aligned to odor onset, with vertical reference
%      lines indicating odor onset (dotted black), reward timing (dashed blue),
%      and airpuff timing (dashed red)
%    - Mean waveform
%    - Autocorrelation function
%    - Log interspike interval (ISI) histogram
%    - Session-wide firing rate time series
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
folder_list = {'D:\vCA1\M593\preprocessed\M593_2025_12_27'};
for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});
    % Loading data
    rng(2024);
    LoadExpKeys;

    % One probe obviously exists
    cfg_master = []; % holds config options for later in this script
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
        S.label = cat(2,S.label,S2.label);
        S.t = cat(2,S.t,S2.t);
        S.usr.recording_depths = cat(2,S.usr.recording_depths, S2.usr.recording_depths);
        S.usr.shank_ids = cat(2,S.usr.shank_ids, S2.usr.shank_ids);
        S.usr.mean_waveforms = cat(2,S.usr.mean_waveforms, S2.usr.mean_waveforms);
        if (isfield(S.usr, 'channel_ids')) & (isfield(S2.usr, 'channel_ids'))
            S.usr.channel_ids = cat(2,S.usr.channel_ids,S2.usr.channel_ids);
        end
        clear S2
    end

    % Making binned version
    cfg_Q = [];
    cfg_Q.dt = 1;
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_size = 5; % gaussian window size in seconds
    cfg_Q.gausswin_sd = 1; % SD for gaussian convolution in seconds
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
    cfg_MUA.tvec = csc.tvec; % timebase to compute MUA on
    MUA = getMUA(cfg_MUA, S);
    MUAz = zscore_tsd(MUA);

    %%  Load events and classify what odors are what
    % Load events (in this case, odor onsets)
    load('all_events.mat')

    % Define odor categories based on ExpKeys
    % ExpKeys should contain fields like:
    % ExpKeys.reward_odors = {'odor1', 'odor2'}
    % ExpKeys.airpuff_odors = {'odor3', 'odor4'}
    % ExpKeys.neutral_odors = {'odor5', 'odor6'}

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

    %% Load unconditioned stimulus events + calculate CS-US latency
    % Extract reward and airpuff events from the events file
    reward_times = evt.t(strcmp(evt.label, 'Reward ON'));
    airpuff_times = evt.t(strcmp(evt.label, 'Airpuff ON'));

    % Extract numeric values from cells
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

    % Calculate mean US latencies from second non-neutral odor onset
    us_window = [1 2]; % US expected in this window after odor
    reward_latencies = [];
    airpuff_latencies = [];

    % Get times for reward and airpuff odors only (non-neutral)
    reward_odor_times = evt.t(reward_indices);
    airpuff_odor_times = evt.t(airpuff_indices);

    % Extract numeric values
    if ~isempty(reward_odor_times)
        reward_odor_times = cat(1, reward_odor_times{:});
    else
        reward_odor_times = [];
    end

    if ~isempty(airpuff_odor_times)
        airpuff_odor_times = cat(1, airpuff_odor_times{:});
    else
        airpuff_odor_times = [];
    end

    % Calculate reward US latencies
    for i = 1:length(reward_times)
        % Find reward odor pairs preceding this US
        preceding_reward_odors = reward_odor_times(reward_odor_times < reward_times(i));

        if length(preceding_reward_odors) >= 2
            % Get the second reward odor in the pair
            second_odor_time = preceding_reward_odors(end);
            latency = reward_times(i) - second_odor_time;

            if latency > us_window(1) && latency < us_window(2)
                reward_latencies = [reward_latencies, latency];
            end
        end
    end

    % Calculate airpuff US latencies
    for i = 1:length(airpuff_times)
        % Find airpuff odor pairs preceding this US
        preceding_airpuff_odors = airpuff_odor_times(airpuff_odor_times < airpuff_times(i));

        if length(preceding_airpuff_odors) >= 2
            % Get the second airpuff odor in the pair
            second_odor_time = preceding_airpuff_odors(end);
            latency = airpuff_times(i) - second_odor_time;

            if latency > us_window(1) && latency < us_window(2)
                airpuff_latencies = [airpuff_latencies, latency];
            end
        end
    end

    % Calculate mean latencies
    if ~isempty(reward_latencies)
        mean_reward_latency = mean(reward_latencies);
    else
        mean_reward_latency = NaN;
    end

    if ~isempty(airpuff_latencies)
        mean_airpuff_latency = mean(airpuff_latencies);
    else
        mean_airpuff_latency = NaN;
    end
    %%  Set up colors and PETH parameters
    % Create a color map with enough colors for all odors
    num_cues = 7; % total number of unique odor labels
    cols = distinguishable_colors(num_cues);

    % the 4th color is always black-ish, so switching that with neutral if needed
    if num_cues >= 4
        cols(4,:) = cols(3,:);
        cols(3,:) = [0 0 0];
    end

    % For odor PETH
    cfg_peth = []; % parameters for PETH
    cfg_peth.window = [-2 2];
    cfg_peth.dt = 0.01;
    cfg_peth.mode = 'interp';


    % Create cue events for all odors
    for iC = 1:num_cues
        cue_evt(iC) = ts(evt.t(iC), evt.label(iC));
        % remove illegal events
        tstart = cue_evt(iC).t{1} + cfg_peth.window(1);
        tend = cue_evt(iC).t{1} + cfg_peth.window(2);
        remove_idx = tstart < MUA.tvec(1) | tend > MUA.tvec(end);
        cue_evt(iC).t{1}(remove_idx) = [];
    end

    %% Create US event objects
    reward_evt = ts({reward_times});
    airpuff_evt = ts({airpuff_times});
    %% Generate cell fingerprints
    for iC = 1:length(S.t)
        this_fig = figure('WindowState', 'maximized');
        this_cell = SelectTS([],S,iC);

        % Calculate firing rate and z-scored firing rate for this cell
        MUA_cell = getMUA(cfg_MUA, this_cell);
        MUAz_cell = zscore_tsd(MUA_cell);

        % Plot odor-PETH (row 1, columns 1-2)
        ax1 = subplot(3,5,[3,4]); % Odor Raw PETH
        hold on
        ax2 = subplot(3,5,[1,2]); % Odor Z-scored PETH
        hold on
        
        % Plot US-PETH (row 2, columns 1-2)
        ax3 = subplot(3,5,[6,7]); % US Raw PETH
        hold on
        ax4 = subplot(3,5,[8,9]); % US Z-scored PETH
        hold on

        leg_labels = {};
        min_ylim = 10;
        max_ylim = -10;

        % Plot all odors
        for iCue = 1:num_cues
            this_cue = restrict(cue_evt(iCue), iv(ExpKeys.sessstart, ExpKeys.sessend));

            if ~isempty(this_cue.t{1})
                this_out = TSDpeth_fast(cfg_peth, MUAz_cell, this_cue);
                this_out_raw = TSDpeth_fast(cfg_peth, MUA_cell, this_cue);
                leg_labels{length(leg_labels)+1} = cue_evt(iCue).label{1};

                % Use a unique color for each odor
                plot(ax1, this_out.tvec, this_out.data, 'Color', cols(iCue,:), 'LineWidth', 1.5);
                min_ylim = min(min_ylim, min(this_out.data)-0.1);
                max_ylim = max(max_ylim, max(this_out.data)+0.1);
                plot(ax2, this_out_raw.tvec, this_out_raw.data, 'Color', cols(iCue,:), 'LineWidth', 1.5);
            end
        end

        % Add US indicator lines for odor PETHs
        xline(ax1, 0, ':', 'Color', [0 0 0], 'LineWidth', 1.5);
        xline(ax2, 0, ':', 'Color', [0 0 0], 'LineWidth', 1.5);

        if ~isnan(mean_airpuff_latency)
            xline(ax1, mean_airpuff_latency, '--', 'Color', [1 0 0], 'LineWidth', 1.5, 'Label', 'Airpuff');
            xline(ax2, mean_airpuff_latency, '--', 'Color', [1 0 0], 'LineWidth', 1.5, 'Label', 'Airpuff');
        end

        if ~isnan(mean_reward_latency)
            xline(ax1, mean_reward_latency, '--', 'Color', [0 0 1], 'LineWidth', 1.5, 'Label', 'Reward');
            xline(ax2, mean_reward_latency, '--', 'Color', [0 0 1], 'LineWidth', 1.5, 'Label', 'Reward');
        end

        % Format raw odor PETH plot
        ax2.Title.String = 'Raw Odor PETH';
        ax2.Title.FontSize = 14;
        ax2.YLabel.String = 'F.R. (Hz)';
        ax2.YLabel.FontSize = 14;
        ax2.XLabel.String = 'Time from odor onset (s)';
        ax2.XLabel.FontSize = 14;

        % Format z-scored odor PETH plot
        ax1.YLim = [min(min_ylim, -0.5), max(max_ylim, 1)];
        ax1.Title.String = 'Z-scored Odor PETH';
        ax1.Title.FontSize = 14;
        ax1.YLabel.String = 'Z-score';
        ax1.YLabel.FontSize = 14;
        ax1.XLabel.String = 'Time from odor onset (s)';
        ax1.XLabel.FontSize = 14;

        % Shared legend for both odor PETHs
        % legend(ax2, leg_labels, 'location', 'best', 'FontSize', 8, 'Box', 'off');

        % Shared legend for both odor PETHs - place to the LEFT
        ax_legend_odor = axes('Position', [0.02, 0.65, 0.08, 0.25]);
        ax_legend_odor.Visible = 'off';
        
        % Create dummy lines for legend
        h_dummy = [];
        for iCue = 1:length(leg_labels)
            h_dummy = [h_dummy, line(NaN, NaN, 'Color', cols(iCue,:), 'LineWidth', 2)];
        end
        legend(ax_legend_odor, h_dummy, leg_labels, 'Location', 'best', 'FontSize', 8, 'Box', 'off');

        % Plot US time-locked PETHs
        leg_labels_us = {};
        min_ylim_us = 10;
        max_ylim_us = -10;

        % Plot reward US PETH
        reward_cue = restrict(reward_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(reward_cue.t{1})
            this_out_us_z = TSDpeth_fast(cfg_peth, MUAz_cell, reward_cue);
            this_out_us_raw = TSDpeth_fast(cfg_peth, MUA_cell, reward_cue);
            leg_labels_us{length(leg_labels_us)+1} = 'Reward';

            plot(ax3, this_out_us_raw.tvec, this_out_us_raw.data, 'Color', [0 0 1], 'LineWidth', 2);
            plot(ax4, this_out_us_z.tvec, this_out_us_z.data, 'Color', [0 0 1], 'LineWidth', 2);

            min_ylim_us = min(min_ylim_us, min(this_out_us_z.data)-0.1);
            max_ylim_us = max(max_ylim_us, max(this_out_us_z.data)+0.1);
        end

        % Plot airpuff US PETH
        airpuff_cue = restrict(airpuff_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));
        if ~isempty(airpuff_cue.t{1})
            this_out_us_z = TSDpeth_fast(cfg_peth, MUAz_cell, airpuff_cue);
            this_out_us_raw = TSDpeth_fast(cfg_peth, MUA_cell, airpuff_cue);
            leg_labels_us{length(leg_labels_us)+1} = 'Airpuff';

            plot(ax3, this_out_us_raw.tvec, this_out_us_raw.data, 'Color', [1 0 0], 'LineWidth', 2);
            plot(ax4, this_out_us_z.tvec, this_out_us_z.data, 'Color', [1 0 0], 'LineWidth', 2);

            min_ylim_us = min(min_ylim_us, min(this_out_us_z.data)-0.1);
            max_ylim_us = max(max_ylim_us, max(this_out_us_z.data)+0.1);
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

        % Shared legend for both US PETHs
        % legend(ax3, leg_labels_us, 'location', 'best', 'FontSize', 10);
        % Shared legend for both US PETHs - place to the LEFT
        ax_legend_us = axes('Position', [0.02, 0.35, 0.08, 0.25]);
        ax_legend_us.Visible = 'off';
        
        % Create dummy lines for legend
        h_dummy_us = [];
        for iUS = 1:length(leg_labels_us)
            if strcmp(leg_labels_us{iUS}, 'Reward')
                h_dummy_us = [h_dummy_us, line(NaN, NaN, 'Color', [0 0 1], 'LineWidth', 2)];
            elseif strcmp(leg_labels_us{iUS}, 'Airpuff')
                h_dummy_us = [h_dummy_us, line(NaN, NaN, 'Color', [1 0 0], 'LineWidth', 2)];
            end
        end
        legend(ax_legend_us, h_dummy_us, leg_labels_us, 'Location', 'best', 'FontSize', 8, 'Box', 'off');

        % Plot waveform (row 1, column 3)
        subplot(3,5,5)
        [~, max_idx] = max((max(this_cell.usr.mean_waveforms{1},[],2)));
        plot(this_cell.usr.mean_waveforms{1}(max_idx,:))
        title('Mean Waveform', 'FontSize', 14)
        xticks([])
        yticks([])

        % Plot autocorrelation (row 2, column 3)
        subplot(3,5,10)
        cfg_xc = [];
        cfg_xc.smooth = 1;
        cfg_xc.max_t = 1;
        [test_xc, test_t] = ccf(cfg_xc,this_cell.t{1},this_cell.t{1});
        plot(test_t, test_xc);
        xticks(0:1/8:1);
        xlim([0 1])
        xlabel('Time (s)', 'FontSize', 14)
        title('Auto-correlation', 'FontSize', 14)

        % Plot log ISI (row 3, column 3)
        subplot(3,5,15)
        test_isi = diff(this_cell.t{1});
        histogram(log10(test_isi), -3.75:0.5:2.75);
        xline(log10(0.002), 'r');
        this_str = sprintf(', %d/%d ISIs less than 2 msec', ...
            sum(test_isi<0.002), length(test_isi));
        title(strcat('LogISI', this_str), 'FontSize', 14)
        ax = gca;
        ax.YAxis.Exponent = 0;
        xlabel('Time (ms)', 'FontSize', 14)

        % Plot firing-rate over time (row 3, columns 1-2)
        subplot(3,5,[11:14])
        hold on
        plot(Q.tvec, Q.data(iC,:), 'black');
        xlabel('Time (s)', 'FontSize', 14)
        ylabel('F.R. (Hz)', 'FontSize', 14)
        xlim([0, ExpKeys.sessend])
        title('Firing Rate Over Time', 'FontSize', 14)

        fname = strcat(ExpKeys.subject, '_', strrep(ExpKeys.date, '-', '_'), ...
            '_', this_cell.label{1}, '_cell_fingerprint_odorUS');
        sgtitle(fname,'Interpreter', 'None');
        print(this_fig, '-dpng',  fname)
        % print(this_fig, '-depsc',  fname)
        close;
    end


    clearvars -except folder_list iFolder
end