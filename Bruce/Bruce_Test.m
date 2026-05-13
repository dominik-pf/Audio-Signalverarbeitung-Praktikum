clear; close;
scriptFolder = fileparts(mfilename('fullpath'));
modelFolder = fullfile(scriptFolder, 'Model');
addpath(genpath(modelFolder));
hubPath = fileparts(scriptFolder);
addpath(genpath(hubPath));

%% model run settings
% Check to see if running under Matlab or Octave
if exist ('OCTAVE_VERSION', 'builtin') ~= 0
  pkg load signal;
  if exist('rms')<1
    rms = @(x) sqrt(mean(x.^2));
  end
end

if exist('parfor','builtin') % check if the Matlab Parallel Computation
                             % Toolbox is installed and use appropriate
                             % function
    generate_neurogram_function = @generate_neurogram_BEZ2018a_parallelized;
    disp('Using parallelized version of neurogram generation function')
else
    generate_neurogram_function = @generate_neurogram_BEZ2018a;
    disp('Using serial version of neurogram generation function')
end

% Set audiogram
ag_fs = [125 250 500 1e3 2e3 4e3 8e3];
ag_dbloss = [0 0 0 0 0 0 0]; % Normal hearing

species = 2; % Human cochlear tuning (Shera et al., 2002)


%% Select input
nr = 0; % select input by changing nr

if nr==0
    % create a custom input with different frequencies
    Fs_stim = 100000;
    t_stim = [0:1/Fs_stim:1]';
    stim = zeros(size(t_stim));

    ix1 = t_stim < 0.3;
    stim(ix1) = stim(ix1) + sin(2*pi*200*t_stim(ix1));
    ix2 = t_stim > 0.2 & t_stim < 0.5;
    stim(ix2) = stim(ix2) + sin(2*pi*4000*t_stim(ix2));
    ix3 = t_stim > 0.4 & t_stim < 0.8;
    stim(ix3) = stim(ix3) + sin(2*pi*1200*t_stim(ix3));
    ix4 = t_stim > 0.8 & t_stim < 1;
    stim(ix4) = stim(ix4) + sin(2*pi*500*t_stim(ix4));

elseif nr==1
    [stim, Fs_stim] = audioread('long_test.wav');
    % only use a segment of the audio
    max_dur = 1.5;   % seconds
    N = min(length(stim), round(max_dur*Fs_stim));
    stim = stim(1:N);
    stim = stim(:);

elseif nr==2
    [stim, Fs_stim] = audioread('defineit.wav');
end

stimdb = 65; % speech level in dB SPL

stim = stim/rms(stim)*20e-6*10^(stimdb/20);

%% Run Bruce model
[neurogram_ft,neurogram_mr,neurogram_Sout,t_ft,t_mr,t_Sout,CFs] = generate_neurogram_function(stim,Fs_stim,species,ag_fs,ag_dbloss);


%%
winlen = 256; % Window length for the spectrogram analyses

%% Meanrate
ng1=figure;
set(ng1,'renderer','painters');

sp1 = subplot(2,1,1);
[s,f,t] = specgram([stim; eps*ones(round(t_mr(end)*Fs_stim)-length(stim),1)],winlen,Fs_stim,winlen,0.25*winlen);
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
xlabel('Time');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;

sp2=subplot(2,1,2);
plot_neurogram(t_mr,CFs,neurogram_mr,sp2);
caxis([0 80])
title('Mean-rate Neurogram')
xlim(xl)

%% S_out
ng3=figure;
set(ng3,'renderer','painters');

sp1 = subplot(2,1,1);
[s,f,t] = specgram([stim; eps*ones(round(t_mr(end)*Fs_stim)-length(stim),1)],winlen,Fs_stim,winlen,0.25*winlen);
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
xlabel('Time');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;

sp2=subplot(2,1,2);
plot_neurogram(t_Sout,CFs,neurogram_Sout*diff(t_Sout(1:2)),sp2);
caxis([0 6])
title('S_{out} Neurogram')
xlim(xl)
