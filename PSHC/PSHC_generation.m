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
idx = 3;                 % 1 = 98 pps, 2 = 200 pps, 3 = 450 pps
envRate = envRates(idx);
k = kVals(idx);

%% Harmonic range covering the passband
M = 1000%1000;%ceil(250 / f0);
N = 1500%1270;%floor(8000 / f0);
harmonics = M:N;

%% Generate PSHCs
% my code
downPSHC = generate_pshc(t, f0, harmonics, k, "down");
upPSHC   = generate_pshc(t, f0, harmonics, k, "up");

% original code
[out_gen_down,t_gen_down]=GenPSHC(f0,M,N,k,dur,fs,'down');
[out_gen_up,t_gen_up]=GenPSHC(f0,M,N,k,dur,fs,'up');

%% Bandpass filter
[b,a] = butter(3, [fLow fHigh]/(fs/2), "bandpass"); % should be 6

downPSHC = filtfilt(b,a,downPSHC);
upPSHC   = filtfilt(b,a,upPSHC);

out_gen_down = filtfilt(b,a,out_gen_down);
out_gen_up = filtfilt(b,a,out_gen_up);

% Apply 20-ms raised-cosine ramps
downPSHC = apply_ramp(downPSHC, fs, rampDur);
upPSHC   = apply_ramp(upPSHC, fs, rampDur);

out_gen_down = apply_ramp(out_gen_down, fs, rampDur);
out_gen_up   = apply_ramp(out_gen_up, fs, rampDur);


%% Plot waveforms and Hilbert envelopes
figure;

subplot(2,1,1)
plot(t, downPSHC); hold on;
plot(t, abs(hilbert(downPSHC)), 'r');
title("Down-PSHC")
xlim([0 0.05])

subplot(2,1,2)
plot(t, upPSHC); hold on;
plot(t, abs(hilbert(upPSHC)), 'r');
title("Up-PSHC")
xlim([0 0.05])

% original
figure;

subplot(2,1,1)
plot(t, out_gen_down); hold on;
plot(t, abs(hilbert(out_gen_down)), 'r');
title("Down-PSHC_orig")
xlim([0 0.05])

subplot(2,1,2)
plot(t, out_gen_up); hold on;
plot(t, abs(hilbert(out_gen_up)), 'r');
title("Up-PSHC_orig")
xlim([0 0.05])

%% Spectrograms
figure;

subplot(1,2,1)
spectrogram(downPSHC, 1024, 900, 2048, fs, "yaxis");
title("Down-PSHC spectrogram")
ylim([1 4])
xlim([0 0.5])
colormap(jet)
colorbar
clim([-100 -20])

subplot(1,2,2)
spectrogram(upPSHC, 1024, 900, 2048, fs, "yaxis");
title("Up-PSHC spectrogram")
ylim([1 4])
xlim([0 0.5])
colormap(jet)
colorbar
clim([-100 -20])

% original
figure;

subplot(1,2,1)
spectrogram(out_gen_down, 1024, 900, 2048, fs, "yaxis");
title("Down-PSHC spectrogram_orig")
ylim([1 4])
xlim([0 0.5])
colormap(jet)
colorbar
clim([-100 -20])

subplot(1,2,2)
spectrogram(out_gen_up, 1024, 900, 2048, fs, "yaxis");
title("Up-PSHC spectrogram_orig")
ylim([1 4])
xlim([0 0.5])
colormap(jet)
colorbar
clim([-100 -20])