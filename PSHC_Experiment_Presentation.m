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
    [downPSHC, upPSHC, fs, t, Frequ_range, envRates] = PSHC_Experiment_Generation_Presentation(1, [250, 500, 1000, 2000, 4000, 8000, 11200]);
end
all_combinations = fieldnames(downPSHC);

%% ::::::::::::::::::::::::::::::::::::::::::
% test = all_combinations{11};
%  :::::::::::::::::::::::::::::::::::::::::: 

%% ::::::::::::::::::::::::::::::::::::::::::
plot_active = true; % plotting on/off
%  ::::::::::::::::::::::::::::::::::::::::::

UP_DOWN_diff = struct();
length(all_combinations)
for combos=1:length(all_combinations)
combos
    test = all_combinations{combos};
    
    freq_range = Frequ_range.(test);
    fLow = freq_range(1);
    fHigh = freq_range(2);

    %% ========================================================================
    %                               downPSHC
    % =========================================================================
    stim_down = downPSHC.(test);
    Fs_stim_down = fs;
    
    stimdb = 65; % speech level in dB SPL
    stim_down = stim_down/rms(stim_down)*20e-6*10^(stimdb/20);
    t_stim_down = t(:);
    stim_down = stim_down(:); % Ensure stim is a column vector
    
    %% Rum Carfac Model
    CF_down = CARFAC_Design(1, Fs_stim_down);
    CF_down = CARFAC_Init(CF_down);
    
    [CF_down, ~, naps_down, ~, ~, ~] = CARFAC_Run(CF_down, stim_down);
    
    
    CFs_down = CF_down.pole_freqs(:);

    if size(naps_down,1) ~= length(t_stim_down)
        naps_down = naps_down.';
    end
    
    cf_down_mask = CFs_down >= 125 & CFs_down <= 20000;
    CFs_down_plot = CFs_down(cf_down_mask);
    naps_down_plot = naps_down(:,cf_down_mask);
    
    %% ========================================================================
    %                               upPSHC
    % =========================================================================
    stim_up = upPSHC.(test);
    Fs_stim_up = fs;
    
    stimdb = 65; % speech level in dB SPL
    stim_up = stim_up/rms(stim_up)*20e-6*10^(stimdb/20);
    t_stim_up = t(:);
    stim_up = stim_up(:); % Ensure stim is a column vector
    
    %% Rum Carfac Model
    CF_up = CARFAC_Design(1, Fs_stim_up);
    CF_up = CARFAC_Init(CF_up);
    
    [CF_up, ~, naps_up, ~, ~, ~] = CARFAC_Run(CF_up, stim_up);
    
    
    CFs_up = CF_up.pole_freqs(:);
    
    if size(naps_up,1) ~= length(t_stim_up)
        naps_up = naps_up.';
    end

    cf_up_mask = CFs_up >= 125 & CFs_up <= 20000; 
    CFs_up_plot = CFs_up(cf_up_mask);
    naps_up_plot = naps_up(:,cf_up_mask);

    %% Cut the signal beginning
    cutTime = 0.1;
    fadeTime = 0.03; % 30 ms smooth transition
    
    t_resp = (0:size(naps_down,1)-1)'/fs;
    
    w_time = ones(size(t_resp));
    w_time(t_resp < cutTime) = 0;
    
    fadeIdx = t_resp >= cutTime & t_resp < cutTime + fadeTime;
    w_time(fadeIdx) = 0.5 - 0.5*cos(pi*(t_resp(fadeIdx)-cutTime)/fadeTime);

    naps_down = naps_down .* w_time;
    naps_up   = naps_up   .* w_time;
  
    %% =========================================================================
    %               Model Analysis
    %  =========================================================================
    tokens = regexp(test, 'Fc(\d+)_env(\d+)', 'tokens', 'once');
    FcField  = ['Fc'  tokens{1}];
    envField = ['env' tokens{2}];

    %% local analysis
    
    % hard cut frequency band
        Fc = str2double(tokens{1});
        CFs = CFs_up;
        [~, ch0] = min(abs(CFs - Fc));
        idx = max(1,ch0-15) : min(length(CFs),ch0+15);

        R_up   = naps_up(:,idx)';
        R_down = naps_down(:,idx)';
        CF_local = CFs(idx);

    % % channel weighting around Fc
    %     Fc = str2double(tokens{1});
    %     CFs = CFs_up(:);
    % 
    %     sigma_oct = 0.55; % bandwidth in octaves, tune this
    %     dist_oct = log2(CFs / Fc);
    % 
    %     w_cf = exp(-0.5 * (dist_oct / sigma_oct).^2);
    %     w_cf = w_cf(:)';
    % 
    %     R_down = naps_down .* w_cf;
    %     R_up   = naps_up   .* w_cf;
    %     R_down = R_down';
    %     R_up = R_up';

    % Structural similarity (SSIM) index for measuring image quality
    [SSIM_score_local, SSIM_diff_map_local] = ssim(R_down, R_up);
    
                                figure
                                imagesc(t_stim_down, CFs_up, SSIM_diff_map_local)
                                axis xy
                                xlabel('Time (s)')
                                ylabel('CF (Hz)')
                                title(sprintf('SSIM difference map: %s', test))
                                colorbar
    
    UP_DOWN_SSIM_score_local.(FcField).(envField) = 1- SSIM_score_local;

    % NSIM 
    A_L = R_down;
    B_L = R_up;
    
    % make sure both are time x CF
    if size(A_L,1) ~= size(B_L,1)
        B_L = B_L.';
    end
    
    [nsim_val_local, nsim_map_local] = nsim_paper(A_L, B_L);
    UP_DOWN_NSIM_score_local.(FcField).(envField) = 1- nsim_val_local;

                                figure
                                imagesc(t_stim_down, CFs_up, nsim_map_local)
                                axis xy
                                xlabel('Time (s)')
                                ylabel('CF (Hz)')
                                title(sprintf('NSIM difference map: %s', test))
                                colorbar

    % %% Analysis between input signals
    %     nwin = 512;
    %     noverlap = 256;
    %     nfft = 2*1024;
    % 
    %     [s_down,f_down,t_spec] = spectrogram(stim_down, nwin, noverlap, nfft, Fs_stim_down, "yaxis");
    %     [s_up,f_up,~] = spectrogram(stim_up, nwin, noverlap, nfft, Fs_stim_up, "yaxis");
    %     s_down = abs(s_down);
    %     s_up = abs(s_up);
    % 
    %     % hard cur frequency band
    %         % mask = (f_down >= CF_local(end)) & (f_down <= CF_local(1));
    %         % s_down = s_down(mask,:);
    %         % f_down = f_down(mask);
    %         % s_up = s_up(mask,:);
    %         % f_up = f_up(mask);
    % 
    %     % frequency weighting around Fc
    %         dist_oct_spec = log2(f_down / Fc);
    %         w_f = exp(-0.5 * (dist_oct_spec / sigma_oct).^2);
    % 
    %         s_down = s_down .* w_f;
    %         s_up   = s_up   .* w_f;
    % 
    %     % Structural similarity (SSIM) index for measuring image quality
    %     [SSIM_score_input, SSIM_diff_map_in] = ssim(s_down, s_up);
    % 
    %                             % figure
    %                             % imagesc(t_stim_down, CFs_up, SSIM_diff_map_in)
    %                             % axis xy
    %                             % xlabel('Time (s)')
    %                             % ylabel('CF (Hz)')
    %                             % title(sprintf('SSIM difference map: %s', test))
    %                             % colorbar
    % 
    %     UP_DOWN_SSIM_score_input.(FcField).(envField) = 1- SSIM_score_input;
    % 
    %     % NSIM 
    %     A_input = s_down;
    %     B_input = s_up;
    % 
    %     % make sure both are time x CF
    %     if size(A_input,1) ~= size(B_input,1)
    %         B_input = B_input.';
    %     end
    % 
    %     [nsim_val_input, nsim_map_input] = nsim_paper(A_input, B_input);
    %     UP_DOWN_NSIM_score_input.(FcField).(envField) = 1- nsim_val_input;
    % 
    %                             % figure
    %                             % imagesc(t_stim_down, CFs_up, nsim_map_input)
    %                             % axis xy
    %                             % xlabel('Time (s)')
    %                             % ylabel('CF (Hz)')
    %                             % title(sprintf('NSIM difference map input: %s', test))
    %                             % colorbar


        % %% Test Alignment
        % maxLag_ms = 10;
        % maxLag = round(maxLag_ms/1000 * fs);
        % 
        % lags = -maxLag:maxLag;
        % diffVals = zeros(size(lags));
        % corrVals = zeros(size(lags));
        % 
        % for k = 1:numel(lags)
        %     lag = lags(k);
        % 
        %     if lag > 0
        %         A = R_down(:,1+lag:end);
        %         B = R_up(:,1:end-lag);
        %     elseif lag < 0
        %         A = R_down(:,1:end+lag);
        %         B = R_up(:,1-lag:end);
        %     else
        %         A = R_down;
        %         B = R_up;
        %     end
        %     [nsim_val, ~] = nsim_paper(A, B);
        %     diffVals(k) = 1 - nsim_val;
        % 
        %     corrVals(k) = corr(A(:), B(:), 'Rows', 'complete');
        % end
        % 
        % [bestDiff, bestIdx] = min(diffVals);
        % bestLag_samples = lags(bestIdx);
        % bestLag_ms = bestLag_samples / fs * 1000;
        % 
        % [bestCorr, bestCorrIdx] = max(corrVals);
        % bestCorrLag_samples = lags(bestCorrIdx);
        % bestCorrLag_ms = bestCorrLag_samples / fs * 1000;
        % 
        % UP_DOWN_NSIM_score_local_min.(FcField).(envField) = bestDiff;
        % UP_DOWN_bestLag_ms.(FcField).(envField) = bestLag_ms;
        % 
        % UP_DOWN_corr_local_max.(FcField).(envField) = bestCorr;
        % UP_DOWN_corr_bestLag_ms.(FcField).(envField) = bestCorrLag_ms;
        
end

toc


%% plot score over envRate NSIM_LOCAL
fcFields = fieldnames(UP_DOWN_NSIM_score_local);

figure; hold on

for f = 1:numel(fcFields)
    fcName = fcFields{f};
    envFields = fieldnames(UP_DOWN_NSIM_score_local.(fcName));

    envRate = zeros(numel(envFields),1);
    diffScore = zeros(numel(envFields),1);

    for i = 1:numel(envFields)
        name = envFields{i};  
        envRate(i) = sscanf(name, 'env%d');
        diffScore(i) = UP_DOWN_NSIM_score_local.(fcName).(name);
    end

    [envRate, idx] = sort(envRate);
    diffScore = diffScore(idx);


    gprMdl = fitrgp(envRate, diffScore);
    envFine = linspace(min(envRate),max(envRate),500)';
    diffFit = predict(gprMdl,envFine);

    UP_DOWN_NSIM_fitting.(fcName).envFine = envFine;
    UP_DOWN_NSIM_fitting.(fcName).diffFit = diffFit;

    % fitobj = fit(envRate, diffScore, 'smoothingspline', 'SmoothingParam',1);
    % envFine = linspace(min(envRate), max(envRate), 500);
    % diffFit = feval(fitobj, envFine);

    plot(envRate, diffScore, '-o', 'LineWidth', 2, 'DisplayName', fcName)
    hold on
    plot(envFine,diffFit,'r','LineWidth',2)

end
legend('Measurements','Smoothing spline')
xlabel('Envelope rate')
ylabel('Differnece score')
title('Differnece score over envelope rate NSIM LOCAL')
legend('Location','best')
grid on
hold off






%% plot score over envRate SSIM_LOCAL
fcFields = fieldnames(UP_DOWN_SSIM_score_local);

figure; hold on

for f = 1:numel(fcFields)
    fcName = fcFields{f};
    envFields = fieldnames(UP_DOWN_SSIM_score_local.(fcName));

    envRate = zeros(numel(envFields),1);
    diffScore = zeros(numel(envFields),1);

    for i = 1:numel(envFields)
        name = envFields{i};  
        envRate(i) = sscanf(name, 'env%d');
        diffScore(i) = UP_DOWN_SSIM_score_local.(fcName).(name);
    end

    [envRate, idx] = sort(envRate);
    diffScore = diffScore(idx);

    
    gprMdl = fitrgp(envRate, diffScore);
    envFine = linspace(min(envRate),max(envRate),500)';
    diffFit = predict(gprMdl,envFine);

    % fitobj = fit(envRate, diffScore, 'smoothingspline');
    % envFine = linspace(min(envRate), max(envRate), 500);
    % diffFit = feval(fitobj, envFine);

    UP_DOWN_SSIM_fitting.(fcName).envFine = envFine;
    UP_DOWN_SSIM_fitting.(fcName).diffFit = diffFit;

    plot(envRate, diffScore, '-o', 'LineWidth', 2, 'DisplayName', fcName)
    hold on
    plot(envFine,diffFit,'r','LineWidth',2)
end

xlabel('Envelope rate')
ylabel('Differnece score')
title('Differnece score over envelope rate SSIM LOCAL')
legend('Location','best')
grid on
hold off

% %% plot score over envRate NSIM_LOCAL_min
% fcFields = fieldnames(UP_DOWN_NSIM_score_local_min);
% 
% figure; hold on
% 
% for f = 1:numel(fcFields)
%     fcName = fcFields{f};
%     envFields = fieldnames(UP_DOWN_NSIM_score_local_min.(fcName));
% 
%     envRate = zeros(numel(envFields),1);
%     diffScore = zeros(numel(envFields),1);
% 
%     for i = 1:numel(envFields)
%         name = envFields{i};  
%         envRate(i) = sscanf(name, 'env%d');
%         diffScore(i) = UP_DOWN_NSIM_score_local_min.(fcName).(name);
%     end
% 
%     [envRate, idx] = sort(envRate);
%     diffScore = diffScore(idx);
% 
%     plot(envRate, diffScore, '-o', 'LineWidth', 2, 'DisplayName', fcName)
% end
% 
% xlabel('Envelope rate')
% ylabel('Differnece score')
% title('Differnece score over envelope rate NSIM LOCAL MINIMUM')
% legend('Location','best')
% grid on
% hold off
% 
% toc
