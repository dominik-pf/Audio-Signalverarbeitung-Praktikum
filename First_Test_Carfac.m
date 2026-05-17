%% Carfac - PSHC Up/Down
%% Add necessary paths
clear; close;
scriptFolder = fileparts(mfilename('fullpath'));
modelFolder = fullfile(scriptFolder, 'Model_matlab');
addpath(genpath(modelFolder));
hubPath = fileparts(scriptFolder);
addpath(genpath(hubPath));
inputPath = fullfile(hubPath,'PSHC');
addpath(genpath(inputPath));

% add Carfac
carfacPath = fullfile(hubPath,'Carfac');
addpath(genpath(carfacPath));
%% Set Parameters for PSHC generation
%% Parameters
fs = 96000;              % Sampling rate
dur = 2;                 % 500 ms
t = 0:1/fs:dur-1/fs;

f0 = 2;                  % Fundamental frequency
envRates = [98 200 450]; % pulses per second
kVals = [7 10 15];       % PSHC orders

fLow = 2000;             % Bandpass lower cutoff
fHigh = 2540;            % Bandpass upper cutoff
rampDur = 0.005;          % 5 ms raised-cosine ramp

%% Choose one condition
idx = 2;                 % 1 = 98 pps, 2 = 200 pps, 3 = 450 pps
envRate = envRates(idx);
k = kVals(idx);

%% Harmonic range covering the passband
M = 1000;%ceil(250 / f0);
N = 1500;%floor(8000 / f0);
harmonics = M:N;

%% Generate PSHCs
% my code
downPSHC = generate_pshc(t, f0, harmonics, k, "down");
upPSHC   = generate_pshc(t, f0, harmonics, k, "up");

%% Bandpass filter
[b,a] = butter(3, [fLow fHigh]/(fs/2), "bandpass"); % should be 6

downPSHC = filtfilt(b,a,downPSHC);
upPSHC   = filtfilt(b,a,upPSHC);

% Apply 20-ms raised-cosine ramps
downPSHC = apply_ramp(downPSHC, fs, rampDur);
upPSHC   = apply_ramp(upPSHC, fs, rampDur);

%% Carfac-v2 Model
%% down
stim = downPSHC;
% stim = upPSHC;
Fs_stim = fs;

stimdb = 65; % speech level in dB SPL
stim = stim/rms(stim)*20e-6*10^(stimdb/20);
t_stim = t(:);
stim = stim(:); % Ensure stim is a column vector

%% Rum Carfac Model
CF = CARFAC_Design(1, Fs_stim);
CF = CARFAC_Init(CF);

[CF, decim_naps, naps, BM, ohc, agc] = CARFAC_Run(CF, stim);


CFs = CF.pole_freqs(:);

if size(BM,1) ~= length(t_stim)
    BM = BM.';
end

if size(naps,1) ~= length(t_stim)
    naps = naps.';
end

cf_mask = CFs >= 125 & CFs <= 20000;

CFs_plot = CFs(cf_mask);
BM_plot = BM(:,cf_mask);
naps_plot = naps(:,cf_mask);

%%
nwin = 512;
noverlap = 256;
nfft = 1024;

%% BM plot
figure;
set(gcf,'renderer','painters');

subplot(2,1,1);
[s,f,t] = spectrogram(stim, nwin, noverlap, nfft, Fs_stim, "yaxis");
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(nwin))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs_plot/1e3) Fs_stim/2e3])])
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;
ylim([1.5 3.25])

sp2 = subplot(2,1,2);
plot_carfac_map(t_stim, CFs_plot, BM_plot, sp2);
title('CARFAC Basilar Membrane Output')
ylabel('CF (kHz)')
caxis([-1 2])
xlim(xl)
ylim(log10([1 4]))

%% NAPS plot
figure;
set(gcf,'renderer','painters');

subplot(2,1,1);
[s,f,t] = spectrogram(stim, nwin, noverlap, nfft, Fs_stim, "yaxis");
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(nwin))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs_plot/1e3) Fs_stim/2e3])])
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;
ylim([1.5 3.25])

sp2 = subplot(2,1,2);
plot_carfac_map(t_stim, CFs_plot, naps_plot, sp2);
title('CARFAC NAPS Output')
ylabel('CF (kHz)')
caxis([-1 2])
xlim(xl)
ylim(log10([1 4]))


%% up
stim = upPSHC;
Fs_stim = fs;

stimdb = 65; % speech level in dB SPL
stim = stim/rms(stim)*20e-6*10^(stimdb/20);
t_stim = t(:);
stim = stim(:); % Ensure stim is a column vector

%% Rum Carfac Model
CF = CARFAC_Design(1, Fs_stim);
CF = CARFAC_Init(CF);

[CF, decim_naps, naps, BM, ohc, agc] = CARFAC_Run(CF, stim);


CFs = CF.pole_freqs(:);

if size(BM,1) ~= length(t_stim)
    BM = BM.';
end

if size(naps,1) ~= length(t_stim)
    naps = naps.';
end

cf_mask = CFs >= 125 & CFs <= 20000;

CFs_plot = CFs(cf_mask);
BM_plot = BM(:,cf_mask);
naps_plot = naps(:,cf_mask);

%%
nwin = 512;
noverlap = 256;
nfft = 1024;

%% BM plot
figure;
set(gcf,'renderer','painters');

subplot(2,1,1);
[s,f,t] = spectrogram(stim, nwin, noverlap, nfft, Fs_stim, "yaxis");
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(nwin))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs_plot/1e3) Fs_stim/2e3])])
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;
ylim([1.5 3.25])

sp2 = subplot(2,1,2);
plot_carfac_map(t_stim, CFs_plot, BM_plot, sp2);
title('CARFAC Basilar Membrane Output')
ylabel('CF (kHz)')
% caxis([-1 2])
xlim(xl)
% ylim(log10([1 4]))

%% NAPS plot
figure;
set(gcf,'renderer','painters');

subplot(2,1,1);
[s,f,t] = spectrogram(stim, nwin, noverlap, nfft, Fs_stim, "yaxis");
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(nwin))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs_plot/1e3) Fs_stim/2e3])])
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;
ylim([1.5 3.25])

sp2 = subplot(2,1,2);
plot_carfac_map(t_stim, CFs_plot, naps_plot, sp2);
title('CARFAC NAPS Output')
ylabel('CF (kHz)')
% caxis([-1 2])
xlim(xl)
% ylim(log10([1 4]))
