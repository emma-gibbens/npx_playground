%% setup
cd('D:\vCA1\M649\preprocessed\M649_2026_04_02');
addpath(genpath('C:\Users\mvdmlab\Documents\egibbens\odor-pixels\analysis\matlab'));
LoadExpKeys;                               % loads ExpKeys structure
cfg.probe = 1;
%% Load LFPs and convert to mvdmlab format

if cfg.probe == 1
    record_depth = ExpKeys.probe1_Depth;
    cfg_md.maxdepth = 1000 * record_depth;
    csc = tsd;
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
else
    csc = tsd;
    load('imec1_clean_lfp.mat');

    csc.tvec = imec1.lfp_tvec';
    csc.data = double(imec1.lfp_traces');
    csc.label = cellstr(string(imec1.channel_ids).strip().extractAfter('ap#')'); % Ugly but works!
    csc.cfg.depths = imec1.depths;
    csc.cfg.shank_ids = imec1.shank_ids;

    Fs = struct;
    Fs.SamplingFrequency = double(imec1.lfp_fs);
    csc.cfg.hdr{1,1} = Fs;
    csc.cfg.depths = (cfg_md.maxdepth - csc.cfg.depths) * cos(deg2rad(ExpKeys.probe2_insertion_roll));
end

% Order LFPs from most dorsal to ventral for sensible plotting in
% MultiRaster
[~, sort_idx] = sort(csc.cfg.depths);
csc.data = csc.data(sort_idx,:);
csc.label = csc.label(sort_idx);
csc.cfg.depths = csc.cfg.depths(sort_idx);
csc.cfg.shank_ids = csc.cfg.shank_ids(sort_idx);
clear sort_idx;

% sorts by shank ID
% [~, sort_idx] = sort(csc.cfg.shank_ids);
% csc.data = csc.data(sort_idx,:);
% csc.label = csc.label(sort_idx);
% csc.cfg.depths = csc.cfg.depths(sort_idx);
% csc.cfg.shank_ids = csc.cfg.shank_ids(sort_idx);
% clear sort_idx;

csc = restrict(csc,0,ExpKeys.sessend);
% csc = restrict(csc,100,1100);

% Build the output struct (oneChannel)
keepIdx = 9;
oneChannel = tsd;

% Keep only the selected rows
oneChannel.tvec = csc.tvec;
oneChannel.data   = csc.data(keepIdx, :);
oneChannel.label  = csc.label(keepIdx);
oneChannel.cfg    = csc.cfg;                 % copy whole cfg
oneChannel.cfg.depths = csc.cfg.depths(keepIdx);     % but trim depths
oneChannel.cfg.shank_ids = csc.cfg.shank_ids(keepIdx); % and shank ids

%% Plot the PSD for the oneChannel
% First decimate
dsf = 2;
oneChannel.data = decimate(oneChannel.data, dsf);
oneChannel.tvec = downsample(oneChannel.tvec, dsf);
oneChannel.cfg.hdr{1}.SamplingFrequency = oneChannel.cfg.hdr{1}.SamplingFrequency./dsf;

Fs = oneChannel.cfg.hdr{1}.SamplingFrequency;

% wSize = 1024;
% [Pxx,F] = periodogram(oneChannel.data,hamming(length(oneChannel.data)),length(oneChannel.data),Fs);
% plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
% xlim([0 150]);

wSize   = 1024*2; % window length (samples)
noverlap = wSize/2;            % 50% overlap between segments
nfft    = 2^nextpow2(wSize);   % FFT length
[Pxx,F] = pwelch(oneChannel.data, hamming(wSize), noverlap, nfft, Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);