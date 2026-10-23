%% Generate PETH for one cell and stimulus onset
%% Input information
cd('D:\vCA1\M593\preprocessed\M593_2025_12_26');

cfg.probe = 1; % change to 1 or 2 depending on which probe the spiking unit is on
cfg.spikeIDs = {'imec0_93'};
cfg.stimulus_label = 'Odor F'; % change to 'Odor A', 'Odor B', 'Reward ON', 'Airpuff ON', etc.
cfg.window = [-1.5 1.5]; % time window around stimulus (s)
cfg.dt = 0.01; % bin size for PETH (s)
cfg.mode = 'interp'; % interpolation mode for PETH
cfg.baseline_window = [-1 0]; % baseline period (s) - relative to stimulus onset
cfg.post_window = [0 1]; % post-stimulus period (s) - relative to stimulus onset

addpath(genpath('C:\Users\mvdmlab\Documents\egibbens\odor-pixels'));
LoadExpKeys; % loads ExpKeys structure

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

% Sort all units according to depth (sensible while plotting)
[~, sort_idx] = sort(S.usr.recording_depths, 'descend');
S.t = S.t(sort_idx);
S.label = S.label(sort_idx);
S.usr.mean_waveforms = S.usr.mean_waveforms(sort_idx);
S.usr.recording_depths = S.usr.recording_depths(sort_idx);
S.usr.shank_ids = S.usr.shank_ids(sort_idx);
clear sort_idx;

%% Load LFPs and convert to mvdmlab format
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

%% Load events
load('all_events.mat')

%% Isolate the spiking data for neurons of interest
S_plot = ts;
if ~isempty(cfg.spikeIDs)
    spikeIdx = find(cellfun(@(x) ismember(x, cfg.spikeIDs), S.label));
    for i = 1:length(spikeIdx)
        S_plot.label{i} = S.label{spikeIdx(i)};
        S_plot.t{1,i} = S.t{1,spikeIdx(i)};
        S_plot.usr.channel_ids{i} = S.usr.channel_ids{spikeIdx(i)};
        S_plot.usr.recording_depths(i) = S.usr.recording_depths(spikeIdx(i));
        S_plot.usr.shank_ids(i) = S.usr.shank_ids(spikeIdx(i));
        S_plot.usr.mean_waveforms{i} = S.usr.mean_waveforms{spikeIdx(i)};
    end
else
    S_plot = S;
end

%% Get stimulus times
stimulus_times = evt.t(strcmp(evt.label, cfg.stimulus_label));
if isempty(stimulus_times) || isempty(stimulus_times{1})
    error(sprintf('Stimulus "%s" not found in events', cfg.stimulus_label));
end
stimulus_times = stimulus_times{1};

%% Set up PETH parameters
cfg_peth = [];
cfg_peth.window = cfg.window;
cfg_peth.dt = cfg.dt;
cfg_peth.mode = cfg.mode;

%% Create stimulus event
stim_evt = ts({stimulus_times});
stim_evt.label = cfg.stimulus_label;

% Remove illegal events (outside recording window)
tstart = stim_evt.t{1} + cfg_peth.window(1);
tend = stim_evt.t{1} + cfg_peth.window(2);
remove_idx = tstart < MUA.tvec(1) | tend > MUA.tvec(end);
stim_evt.t{1}(remove_idx) = [];

%% Generate PETH for each cell in S_plot
for iCell = 1:length(S_plot.t)

    % Select this cell
    this_cell = SelectTS([], S_plot, iCell);

    % Calculate firing rate and z-scored firing rate for this cell
    MUA_cell = getMUA(cfg_MUA, this_cell);
    MUAz_cell = zscore_tsd(MUA_cell);

    % Restrict stimulus events to session window
    this_stim = restrict(stim_evt, iv(ExpKeys.sessstart, ExpKeys.sessend));

    if isempty(this_stim.t{1})
        warning(sprintf('No valid stimulus events for cell %s', this_cell.label{1}));
        continue;
    end

    % Calculate PETHs
    peth_z = TSDpeth_fast(cfg_peth, MUAz_cell, this_stim);
    peth_raw = TSDpeth_fast(cfg_peth, MUA_cell, this_stim);

    %% Calculate baseline and post-stimulus averages
    % For raw PETH
    baseline_idx = peth_raw.tvec >= cfg.baseline_window(1) & peth_raw.tvec < cfg.baseline_window(2);
    post_idx = peth_raw.tvec >= cfg.post_window(1) & peth_raw.tvec <= cfg.post_window(2);

    baseline_fr = mean(peth_raw.data(baseline_idx));
    post_fr = mean(peth_raw.data(post_idx));

    % For z-scored PETH
    baseline_z = mean(peth_z.data(baseline_idx));
    post_z = mean(peth_z.data(post_idx));

    %% Create figure
    fig = figure('WindowState', 'maximized');

    % Plot Z-scored PETH (top)
    ax1 = subplot(2, 1, 1);
    hold on
    plot(ax1, peth_z.tvec, peth_z.data, 'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'Firing Rate');
    xline(ax1, 0, ':', 'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'Stimulus onset');

    % Add baseline and post-stimulus reference lines for z-scored PETH
    plot(ax1, [cfg.window(1), 0], [baseline_z, baseline_z], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', sprintf('Baseline: %.2f', baseline_z));
    plot(ax1, [0, cfg.window(2)], [post_z, post_z], '--', 'Color', [1 0.5 0], 'LineWidth', 1.5, 'DisplayName', sprintf('Post-stim: %.2f', post_z));

    ax1.Title.String = sprintf('Z-scored PETH - %s', cfg.stimulus_label);
    ax1.Title.FontSize = 14;
    ax1.YLabel.String = 'Z-score';
    ax1.YLabel.FontSize = 12;
    ax1.XLabel.String = '';
    grid on;
    legend(ax1, 'Location', 'bestoutside', 'FontSize', 10);

    % Plot Raw PETH (bottom)
    ax2 = subplot(2, 1, 2);
    hold on
    plot(ax2, peth_raw.tvec, peth_raw.data, 'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'Firing Rate');
    xline(ax2, 0, ':', 'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'Stimulus onset');

    % Add baseline and post-stimulus reference lines for raw PETH
    plot(ax2, [cfg.window(1), 0], [baseline_fr, baseline_fr], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', sprintf('Baseline: %.2f Hz', baseline_fr));
    plot(ax2, [0, cfg.window(2)], [post_fr, post_fr], '--', 'Color', [1 0.5 0], 'LineWidth', 1.5, 'DisplayName', sprintf('Post-stim: %.2f Hz', post_fr));

    ax2.Title.String = sprintf('Raw PETH - %s', cfg.stimulus_label);
    ax2.Title.FontSize = 14;
    ax2.YLabel.String = 'Firing Rate (Hz)';
    ax2.YLabel.FontSize = 12;
    ax2.XLabel.String = 'Time from stimulus onset (s)';
    ax2.XLabel.FontSize = 12;
    grid on;
    legend(ax2, 'Location', 'bestoutside', 'FontSize', 10);

    % Title
    fname = strcat(ExpKeys.subject, '_', strrep(ExpKeys.date, '-', '_'), ...
        '_', this_cell.label{1}, '_', cfg.stimulus_label,'_', 'PETH');
    sgtitle(fname, 'Interpreter', 'None', 'FontSize', 14);

    % Save figure
    print(fig, '-dpng', fname);
    close(fig);

    fprintf('Generated PETH for cell %s (Baseline FR: %.2f Hz, Post-stim FR: %.2f Hz)\n', this_cell.label{1}, baseline_fr, post_fr);
end

clearvars -except cfg