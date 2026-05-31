% clear; close all; clc;

%% Parameters
fs = 96000;              % Sampling rate
dur = 2;               % 500 ms
t = 0:1/fs:dur-1/fs;

f0 = 2;                  % Fundamental frequency
envRates = [98 200 450]; % pulses per second
kVals = [7 10 15];       % PSHC orders

fLow = 2000;             % Bandpass lower cutoff
fHigh = 2540;            % Bandpass upper cutoff

rampDur = 0.02;          % 20 ms raised-cosine ramp

%% Choose one condition
idx = 2;                 % 1 = 98 pps, 2 = 200 pps, 3 = 450 pps
envRate = envRates(idx);
k = kVals(idx);

%% Harmonic range covering the passband
M = ceil(100 / f0);
N = floor(20000 / f0);
harmonics = M:N;

%% Generate PSHCs
downPSHC = generate_pshc(t, f0, harmonics, k, "down");
upPSHC   = generate_pshc(t, f0, harmonics, k, "up");

%% Bandpass filter
[b,a] = butter(3, [fLow fHigh]/(fs/2), "bandpass"); % should be 6

downPSHC = filtfilt(b,a,downPSHC);
upPSHC   = filtfilt(b,a,upPSHC);

%% Apply 20-ms raised-cosine ramps
downPSHC = apply_ramp(downPSHC, fs, rampDur);
upPSHC   = apply_ramp(upPSHC, fs, rampDur);

%% listen to the audio
downListen = downPSHC / max(abs(downPSHC));
upListen   = upPSHC   / max(abs(upPSHC));

downListen = 0.5 * downListen;   % 50% volume
upListen   = 0.5 * upListen;

sound(downListen, fs)
pause(3)
sound(upListen, fs)

%% Plot waveforms and Hilbert envelopes
figure
env = abs(hilbert(downPSHC));

plot(t, downPSHC, 'b'); hold on;
plot(t, env, 'r', 'LineWidth', 1.2);
plot(t, -env, 'r', 'LineWidth', 1.2);
xlim([0.045 0.055])
ylim([-0.2 0.2])


figure;

subplot(2,1,1)
plot(t, downPSHC); hold on;
% plot(t, abs(hilbert(downPSHC)), 'r');
title("Down-PSHC")
hold off
% xlim([0 0.05])

subplot(2,1,2)
plot(t, upPSHC); hold on;
% plot(t, abs(hilbert(upPSHC)), 'r');
title("Up-PSHC")
hold off
% xlim([0 0.05])


%% Spectrograms
figure;

subplot(2,1,1)
spectrogram(downPSHC, 1024, 900, 2048, fs, "yaxis");
title("Down-PSHC spectrogram")
ylim([1.5 3])

subplot(2,1,2)
spectrogram(upPSHC, 1024, 900, 2048, fs, "yaxis");
title("Up-PSHC spectrogram")
ylim([1.5 3])
