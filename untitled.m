%% Carfac Envelop Sensitivity
%% Add necessary paths
clear;
% clearvars -except downPSHC upPSHC fs t Frequ_range; 
close all;
scriptFolder = fileparts(mfilename('fullpath'));
modelFolder = fullfile(scriptFolder, 'Model_matlab');
addpath(genpath(modelFolder));
hubPath = fileparts(scriptFolder);
addpath(genpath(hubPath));
inputPath = fullfile(hubPath,'PSHC');
addpath(genpath(inputPath));

tic
% add Carfac
carfacPath = fullfile(hubPath,'Carfac');
addpath(genpath(carfacPath));
 
%% PSHC Generation
if ~exist('downPSHC','var') || isempty(downPSHC)
    [downPSHC, upPSHC, fs, t, Frequ_range, envRates] = PSHC_Experiment_Generation_Presentation(1, 1000);
end
all_combinations = fieldnames(downPSHC);

%% Plot all signals
tests = fieldnames(upPSHC);

for k = 1:length(tests)

    test = tests{k};

    sig_up = upPSHC.(test);
    sig_down = downPSHC.(test);

    nwin = round(0.01 * fs);
    noverlap = round(0.75*nwin);
    nfft = 2^nextpow2(4*nwin);

    figure
    spectrogram(sig_up, hanning(nwin), noverlap, nfft, fs, 'yaxis');
    title('UP PSHC Spectrogram')
    ylim([0 12])
    colorbar

    figure
    spectrogram(sig_down, hanning(nwin), noverlap, nfft, fs, 'yaxis');
    title('DOWN PSHC Spectrogram')
    ylim([0 12])
    colorbar

end