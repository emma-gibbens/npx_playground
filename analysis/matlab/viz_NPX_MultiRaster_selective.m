
%% Paste directory that contains spike and LFP .mat files
cd('D:\Task3-vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting');
LoadExpKeys;
cfg_md = [];
record_depth = ExpKeys.probe1_Depth; %Change depth from this recording. Need to update to accept ExpKeys input.
cfg_md.maxdepth = 1000 * record_depth;

%% Load spikes
% load('clean_units_imec0.mat')
% S = ts;
% S.t  = cellfun(@(x) x', spike_train, 'UniformOutput', false); % This is important to convert the spiketrain into column format
% S.label = cellstr(string(unit_ids).strip())'; % Ugly but works!
% S.usr = [];
% S.usr.recording_depths = depths; % Not sure if this is the best place to save this
% S.usr.mean_waveforms = mean_waveforms; % Not sure if this is the best place to save this
% S.usr.recording_depths = (cfg_md.maxdepth - S.usr.recording_depths) * cos(deg2rad(ExpKeys.probe1_insertion_roll));
% S.usr.shank_ids = shank_ids;

load("selective_permissive_imec0_odor_only.mat");
S = S_out;

% Sort all units according to depth (Sensible while plotting)
[~, sort_idx] = sort(S.usr.odor); % This is because Multi-Raster plots spikes from bottom to top
S.t = S.t(sort_idx);
S.label = S.label(sort_idx);
S.usr.mean_waveforms = S.usr.mean_waveforms(sort_idx);
S.usr.recording_depths = S.usr.recording_depths(sort_idx);
S.usr.shank_ids = S.usr.shank_ids(sort_idx);
clear sort_idx;


%% Load LFPs and convert to mvdmlab format
csc = tsd;

%set for imec0
load('imec0_clean_lfp.mat');

csc.tvec = imec0.lfp_tvec';
csc.data = double(imec0.lfp_traces');
csc.label = cellstr(string(imec0.channel_ids).strip().extractAfter('ap#')'); % Ugly but works!
csc.cfg.depths = imec0.depths;
csc.cfg.shank_ids = imec0.shank_ids;

Fs = struct;
Fs.SamplingFrequency = double(imec0.lfp_fs);
csc.cfg.hdr{1,1} = Fs;
csc.cfg.depths = (cfg_md.maxdepth - csc.cfg.depths) * cos(deg2rad(ExpKeys.probe1_insertion_roll));

% set imec1
% load('imec1_clean_lfp.mat');
% 
% csc.tvec = imec1.lfp_tvec';
% csc.data = double(imec1.lfp_traces');
% csc.label = cellstr(string(imec1.channel_ids).strip().extractAfter('ap#')'); % Ugly but works!
% csc.cfg.depths = imec1.depths;
% csc.cfg.shank_ids = imec1.shank_ids;
% 
% Fs = struct;
% Fs.SamplingFrequency = double(imec1.lfp_fs);
% csc.cfg.hdr{1,1} = Fs;
% csc.cfg.depths = (cfg_md.maxdepth - csc.cfg.depths) * cos(deg2rad(ExpKeys.probe2_insertion_roll));

% Order LFPs from most dorsal to ventral for sensible plotting in
% MultiRaster
[~, sort_idx] = sort(csc.cfg.depths);
csc.data = csc.data(sort_idx,:);
csc.label = csc.label(sort_idx);
csc.cfg.shank_ids = csc.cfg.shank_ids(sort_idx);
csc.cfg.depths = csc.cfg.depths(sort_idx);
clear sort_idx;

[~, sort_idx] = sort(csc.cfg.shank_ids);
csc.data = csc.data(sort_idx,:);
csc.label = csc.label(sort_idx);
csc.cfg.shank_ids = csc.cfg.shank_ids(sort_idx);
csc.cfg.depths = csc.cfg.depths(sort_idx);
clear sort_idx;


%% Prepare LFPs for plotting
cfg_mr = [];

chan_list = 1:1:size(csc.data, 1); 
% chan_list = [25,26];
% chan_list = [28,29,30,31,32,33];
for iCh = 1:length(chan_list)
    this_lfp = csc;
    this_lfp.data = this_lfp.data(chan_list(iCh), :);
    this_lfp.label = this_lfp.label(chan_list(iCh));
    this_lfp.cfg.depths = this_lfp.cfg.depths(chan_list(iCh));
    cfg_mr.lfp(iCh) = this_lfp;
end

%% Prepare events for plotting
load('all_events.mat')
% evt.t = evt.t(1:8);
% evt.label = evt.label(1:8);
% cfg_mr.evt = [];
% if isempty(evt.t{8})
%     evt.t(8) = [];
%     evt.label(8) = [];
% end
% cfg_mr.evt = evt;

evt.t = evt.t(7);
evt.label = evt.label(7);
% cfg_mr.evt = [];
% if isempty(evt.t{8})
%     evt.t(8) = [];
%     evt.label(8) = [];
% end
cfg_mr.evt = evt;


%% Plot MultiRaster
cfg_mr.lfpHeight = 20;
cfg_mr.lfpSpacing = 20;
cfg_mr.lfpMax = 30;
cfg_mr.lfpWidth = 0.25;
cfg_mr.lfpColor = [0,0,0;0,0,0];
% cfg_mr.spkColor = 'w';

% Use this if you want to color the spikes based on their proximity to the
% LFP channels being plotted

% all_colors = linspecer(length(cfg_mr.lfp)); 
% spk_colors = [];
% for iC = 1:length(S.t)
%     [~, min_idx] = min(abs(S.usr.recording_depths(iC) - csc.cfg.depths(chan_list))/...
%         cos(deg2rad(ExpKeys.probe1_insertion_roll))); % being pedantic about the calculation
%     spk_colors = [spk_colors; all_colors(min_idx,:)];
% end
% cfg_mr.spkColor = spk_colors;

% Use this if you want to color only those spikes that are within some
% range of proximity from a given list of channels (say 300 microns)

% all_colors = linspecer(length(cfg_mr.lfp)); 
% spk_colors = [];
% dist_thresh = 300;
% for iC = 1:length(S.t)
%     [min_dist, min_idx] = min(abs(S.usr.recording_depths(iC) - csc.cfg.depths(chan_list))/...
%         cos(deg2rad(ExpKeys.probe1_insertion_roll))); % being pedantic about the calculation
%     if min_dist <= 300
%         spk_colors = [spk_colors; all_colors(min_idx,:)];
%     else
%         spk_colors = [spk_colors; 1,1,1];
%     end
% end
% cfg_mr.spkColor = spk_colors;

A_color = [0,0,1];
B_color = [1,0,0];
C_color = [0,1,0];
D_color = [1,0.5,0.8];
E_color = [0.9290,0.6940,0.1250];
F_color = [0.4940,0.1840,0.5560];

spk_colors = [];

for iC = 1:length(S.t)
    if contains(S.usr.odor(iC),'A') == 1
        spk_colors = [spk_colors;A_color];

    elseif contains(S.usr.odor(iC),'B') == 1
        spk_colors = [spk_colors;B_color];

    elseif contains(S.usr.odor(iC),'C') == 1
        spk_colors = [spk_colors;C_color];

    elseif contains(S.usr.odor(iC),'D') == 1
        spk_colors = [spk_colors;D_color];

    elseif contains(S.usr.odor(iC),'E') == 1
        spk_colors = [spk_colors;E_color];

    elseif contains(S.usr.odor(iC),'F') == 1
        spk_colors = [spk_colors;F_color];

    else sprintf('Unrecognized odor label')
    end

end
cfg_mr.spkColor = spk_colors;


MultiRaster_NPX(cfg_mr, S)
%set(gca, 'Color', [0 0 0]);
% set(gca, 'XColor', [1 1 1]);
% set(gca, 'YColor', [1 1 1]);
%set(gcf, 'Color', [0 0 0]);
set(gcf, 'Renderer','painters')
set(gca, 'YTick', [])
ylabel('')