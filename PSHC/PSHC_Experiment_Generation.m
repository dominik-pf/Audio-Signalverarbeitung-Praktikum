function [downPSHC, upPSHC, fs, t, Frequ_range] = PSHC_Experiment_Generation(envRates, Fc)
%% Parameters
% Input
%   envRates:  select a vector of Envelope Rates that you want to generate

fs = 96000;              % Sampling rate
dur = 2;                 % 500 ms
t = 0:1/fs:dur-1/fs;
rampDur = 0.005;          % 5 ms raised-cosine ramp

f0 = 2;                  % Fundamental frequency
% Fc = [250 500 1000 2000 4000 8000 11200];
% Fc = 8000;

downPSHC = struct();
upPSHC = struct();
Frequ_range = struct();

for idx=1:length(Fc)
    ERBN = (24.7*(4.37*Fc(idx)/1000+1));
    
    fLow = max(0.1,Fc(idx) - ERBN);         % Bandpass lower cutoff
    fHigh = Fc(idx) + ERBN;                 % Bandpass upper cutoff

    %% Harmonic range covering the passband
    M = ceil(100 / f0);
    N = floor(20000 / f0);
    harmonics = M:N;

    for n=1:length(envRates)
        envRate = envRates(n);
        k = round(sqrt(envRate/f0));

        %% Generate PSHCs
        % my code
        downPSHC_it = generate_pshc(t, f0, harmonics, k, "down");
        upPSHC_it   = generate_pshc(t, f0, harmonics, k, "up");
        
        %% Bandpass filter
        [b,a] = butter(3, [fLow fHigh]/(fs/2), "bandpass"); % 6 as order=n*2 for bandpass
        
        downPSHC_it = filtfilt(b,a,downPSHC_it);
        upPSHC_it   = filtfilt(b,a,upPSHC_it);
        
        % Apply 20-ms raised-cosine ramps
        downPSHC_it = apply_ramp(downPSHC_it, fs, rampDur);
        upPSHC_it   = apply_ramp(upPSHC_it, fs, rampDur);

        name = sprintf('Fc%d_env%d', Fc(idx), envRate);
        Frequ_range.(name) = [fLow fHigh]; 
        downPSHC.(name) = downPSHC_it;  % Store down PSHC
        upPSHC.(name) = upPSHC_it;      % Store up PSHC

        clear downPSHC_it upPSHC_it envRate k b a
    end
end