function [downPSHC, upPSHC, fs, t] = PSHC_Experiment_Generation(envRates)
%% Parameters
% Input
%   envRates:  select a vector of Envelope Rates that you want to generate

fs = 96000;              % Sampling rate
dur = 2;                 % 500 ms
t = 0:1/fs:dur-1/fs;
rampDur = 0.005;          % 5 ms raised-cosine ramp

f0 = 2;                  % Fundamental frequency
Fc = [0.25 0.5 1 2 4 8 11.2];

downPSHC = struct();
upPSHC = struct();

for idx=1:length(Fc)
    for n=1:length(envRates)
        ERBN = 24.7*(4.37*Fc(idx)+1);
        
        fLow = max(0.1,Fc(idx) - ERBN);             % Bandpass lower cutoff
        fHigh = Fc(idx) + ERBN;            % Bandpass upper cutoff
        
        envRate = envRates(n);
        k = sqrt(envRate/f0);
        
        %% Harmonic range covering the passband
        M = ceil(fLow / f0);
        N = floor(fHigh / f0);
        harmonics = M:N;
    
        %% Generate PSHCs
        % my code
        downPSHC_it = generate_pshc(t, f0, harmonics, k, "down");
        upPSHC_it   = generate_pshc(t, f0, harmonics, k, "up");
        
        %% Bandpass filter
        [b,a] = butter(3, [fLow fHigh]/(fs/2), "bandpass"); % should be 6
        
        downPSHC_it = filtfilt(b,a,downPSHC_it);
        upPSHC_it   = filtfilt(b,a,upPSHC_it);
        
        % Apply 20-ms raised-cosine ramps
        downPSHC_it = apply_ramp(downPSHC_it, fs, rampDur);
        upPSHC_it   = apply_ramp(upPSHC_it, fs, rampDur);

        name = sprintf('Fc%d_env%d', Fc(idx)*1e3, envRate);
        downPSHC.(name) = downPSHC_it;  % Store down PSHC
        upPSHC.(name) = upPSHC_it;      % Store up PSHC
    end
end