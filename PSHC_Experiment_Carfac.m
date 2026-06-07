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
    [downPSHC, upPSHC, fs, t, Frequ_range, envRates] = PSHC_Experiment_Generation(1, [250, 500, 1000, 2000, 4000, 8000]);
end
all_combinations = fieldnames(downPSHC);

%% ::::::::::::::::::::::::::::::::::::::::::
% test = all_combinations{11};
%  ::::::::::::::::::::::::::::::::::::::::::

%% ::::::::::::::::::::::::::::::::::::::::::
plot_active = true; % plotting on/off
%  ::::::::::::::::::::::::::::::::::::::::::

UP_DOWN_diff = struct();
for combos=1:length(all_combinations)
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
    
    [CF_down, ~, naps_down, BM_down, ~, ~] = CARFAC_Run(CF_down, stim_down);
    
    
    CFs_down = CF_down.pole_freqs(:);
    
    if size(BM_down,1) ~= length(t_stim_down)
        BM_down = BM_down.';
    end
    
    if size(naps_down,1) ~= length(t_stim_down)
        naps_down = naps_down.';
    end
    
    cf_down_mask = CFs_down >= 125 & CFs_down <= 20000;
    CFs_down_plot = CFs_down(cf_down_mask);
    naps_down_plot = naps_down(:,cf_down_mask);

    % Cut the signal beginning
    cutTime = 0.1;                 % seconds
    cutSamples = round(cutTime * Fs_stim_down);
    
    naps_down = naps_down(cutSamples+1:end, :);
    
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
    
    [CF_up, ~, naps_up, BM_up, ~, ~] = CARFAC_Run(CF_up, stim_up);
    
    
    CFs_up = CF_up.pole_freqs(:);
    
    if size(BM_up,1) ~= length(t_stim_up)
        BM_up = BM_up.';
    end
    
    if size(naps_up,1) ~= length(t_stim_up)
        naps_up = naps_up.';
    end

    cf_up_mask = CFs_up >= 125 & CFs_up <= 20000; 
    CFs_up_plot = CFs_up(cf_up_mask);
    naps_up_plot = naps_up(:,cf_up_mask);

    % Cut the signal beginning
    cutTime = 0.1;                 % seconds
    cutSamples = round(cutTime * Fs_stim_down);
    
    naps_up = naps_up(cutSamples+1:end, :);
  
    %% =========================================================================
    %               Model Analysis
    %  =========================================================================
    tokens = regexp(test, 'Fc(\d+)_env(\d+)', 'tokens', 'once');
    FcField  = ['Fc'  tokens{1}];
    envField = ['env' tokens{2}];

    %% local analysis
    Fc = str2double(tokens{1});
    CFs = CFs_up;
    [~, ch0] = min(abs(CFs - Fc));
    idx = max(1,ch0-15) : min(length(CFs),ch0+15);

    R_up   = naps_up(:,idx)';
    R_down = naps_down(:,idx)';
    CF_local = CFs(idx);

    %% Structural similarity (SSIM) index for measuring image quality
    SSIM_score = ssim(naps_down, naps_up);
    
    UP_DOWN_SSIM_score.(FcField).(envField) = 1- SSIM_score;

    %% NSIM 
    A = naps_down;
    B = naps_up;
    
    % make sure both are time x CF
    if size(A,1) ~= size(B,1)
        B = B.';
    end
    
    [nsim_val, ~] = nsim_paper(A, B);
    UP_DOWN_NSIM_score.(FcField).(envField) = 1- nsim_val;


                                            %% Compare only local CFs (around the Fc of PSHC) with "normalized Euclidean distance"
                                            D = norm(R_up(:)-R_down(:)) / norm(R_up(:)+R_down(:));

                                            Channel_diff.(FcField).(envField) = D;

                                            %% Compare only local CFs (around the Fc of PSHC) with correlation
                                            correl = corr(R_up(:), R_down(:));
                                            Channel_correlationDiff.(FcField).(envField) = 1-correl;

    %% local analysis
            % Structural similarity (SSIM) index for measuring image quality
            SSIM_score_local = ssim(R_down, R_up);
            
            UP_DOWN_SSIM_score_local.(FcField).(envField) = 1- SSIM_score_local;
        
            % NSIM 
            A_L = R_down;
            B_L = R_up;
            
            % make sure both are time x CF
            if size(A_L,1) ~= size(B_L,1)
                B_L = B_L.';
            end
            
            [nsim_val_local, ~] = nsim_paper(A_L, B_L);
            UP_DOWN_NSIM_score_local.(FcField).(envField) = 1- nsim_val_local;

        %% Analysis between input signals
            nwin = 512;
            noverlap = 256;
            nfft = 2*1024;
    
            [s_down,f_down,t_spec] = spectrogram(stim_down, nwin, noverlap, nfft, Fs_stim_down, "yaxis");
            [s_up,f_up,~] = spectrogram(stim_up, nwin, noverlap, nfft, Fs_stim_up, "yaxis");
            s_down = abs(s_down);
            s_up = abs(s_up);

            mask = (f_down >= CF_local(end)) & (f_down <= CF_local(1));
            s_down = s_down(mask,:);
            f_down = f_down(mask);
            s_up = s_up(mask,:);
            f_up = f_up(mask);
                        
            % Structural similarity (SSIM) index for measuring image quality
            SSIM_score_input = ssim(s_down, s_up);
     
            UP_DOWN_SSIM_score_input.(FcField).(envField) = 1- SSIM_score_input;
        
            % NSIM 
            A_input = s_down;
            B_input = s_up;
            
            % make sure both are time x CF
            if size(A_input,1) ~= size(B_input,1)
                B_input = B_input.';
            end
            
            [nsim_val_input, ~] = nsim_paper(A_input, B_input);
            UP_DOWN_NSIM_score_input.(FcField).(envField) = 1- nsim_val_input;

        % %% Analysis between input and naps down
        %     % Structural similarity (SSIM) index for measuring image quality
        %     SSIM_score_input_out_down = ssim(s_down, naps_down);
        % 
        %     UP_DOWN_SSIM_score_input_out_down.(FcField).(envField) = 1- SSIM_score_input_out_down;
        % 
        %     % NSIM 
        %     A_input_out_down = s_down;
        %     B_input_out_down = naps_down;
        % 
        %     % make sure both are time x CF
        %     if size(A_input_out_down,1) ~= size(B_input_out_down,1)
        %         B_input_out_down = B_input_out_down.';
        %     end
        % 
        %     [nsim_val_input_out_down, ~] = nsim_paper(A_input_out_down, B_input_out_down);
        %     UP_DOWN_NSIM_score_input_out_down.(FcField).(envField) = 1- nsim_val_input_out_down;
        % 
        % %% Analysis between input and naps up
        %     % Structural similarity (SSIM) index for measuring image quality
        %     SSIM_score_input_out_up = ssim(s_up, naps_up);
        % 
        %     UP_DOWN_SSIM_score_input_out_up.(FcField).(envField) = 1- SSIM_score_input_out_up;
        % 
        %     % NSIM 
        %     A_input_out_up = s_up;
        %     B_input_out_up = naps_up;
        % 
        %     % make sure both are time x CF
        %     if size(A_input_out_up,1) ~= size(B_input_out_up,1)
        %         B_input_out_up = B_input_out_up.';
        %     end
        % 
        %     [nsim_val_input_out_up, ~] = nsim_paper(A_input_out_up, B_input_out_up);
        %     UP_DOWN_NSIM_score_input_out_up.(FcField).(envField) = 1- nsim_val_input_out_up;

        figure
        set(gcf,'renderer','painters');
        subplot(1,2,1);
        imagesc(t_spec,f_down/1e3,20*log10(abs(s_down)/sum(hanning(nwin))*sqrt(2)/20e-6));
        axis xy; axis tight;
        hcb = colorbar;
        set(get(hcb,'ylabel'),'string','SPL')
        caxis([stimdb-80 stimdb])
        xlabel('Time (s)');
        ylabel('Frequency (kHz)');
        title('Spectrogram - PSHC_down - %s',test)
        xl = xlim;

        sp2 = subplot(1,2,2);
        plot_carfac_map(t_stim_down, CF_local, R_down, sp2);
        title('CARFAC NAPS Output - PSHC_down')
        ylabel('CF (kHz)')
        caxis([-1 2])
        xlim(xl)

                                % %%
                                % nwin = 512;
                                % noverlap = 256;
                                % nfft = 2*1024;
                                % 
                                % % NAPS plot down
                                % figure;
                                % set(gcf,'renderer','painters');
                                % 
                                % subplot(2,2,1);
                                % [s,f,t_spec] = spectrogram(stim_down, nwin, noverlap, nfft, Fs_stim_down, "yaxis");
                                % imagesc(t_spec,f/1e3,20*log10(abs(s)/sum(hanning(nwin))*sqrt(2)/20e-6));
                                % axis xy; axis tight;
                                % hcb = colorbar;
                                % set(get(hcb,'ylabel'),'string','SPL')
                                % caxis([stimdb-80 stimdb])
                                % ylim([0 min([max(CFs_down_plot/1e3) Fs_stim_down/2e3])])
                                % xlabel('Time (s)');
                                % ylabel('Frequency (kHz)');
                                % title('Spectrogram - PSHC_down - %s',test)
                                % xl = xlim;
                                % % ylim([fLow fHigh]*1e-3)
                                % 
                                % sp2 = subplot(2,2,3);
                                % plot_carfac_map(t_stim_down, CFs_down_plot, naps_down_plot, sp2);
                                % title('CARFAC NAPS Output - PSHC_down')
                                % ylabel('CF (kHz)')
                                % caxis([-1 2])
                                % xlim(xl)
                                % % ylim(log10([fLow fHigh]*1e-3))
                                % 
                                % % NAPS plot
                                % 
                                % subplot(2,2,2);
                                % [s,f,t_spec] = spectrogram(stim_up, nwin, noverlap, nfft, Fs_stim_up, "yaxis");
                                % imagesc(t_spec,f/1e3,20*log10(abs(s)/sum(hanning(nwin))*sqrt(2)/20e-6));
                                % axis xy; axis tight;
                                % hcb = colorbar;
                                % set(get(hcb,'ylabel'),'string','SPL')
                                % caxis([stimdb-80 stimdb])
                                % ylim([0 min([max(CFs_up_plot/1e3) Fs_stim_up/2e3])])
                                % xlabel('Time (s)');
                                % ylabel('Frequency (kHz)');
                                % title('Spectrogram - PSHC_up')
                                % xl = xlim;
                                % % ylim([fLow fHigh]*1e-3)
                                % 
                                % sp2 = subplot(2,2,4);
                                % plot_carfac_map(t_stim_up, CFs_up_plot, naps_up_plot, sp2);
                                % title('CARFAC NAPS Output - PSHC_up')
                                % ylabel('CF (kHz)')
                                % caxis([-1 2])
                                % xlim(xl)
                                % % ylim(log10([fLow fHigh]*1e-3))

end


%% plot score over envRate NSIM
fcFields = fieldnames(UP_DOWN_NSIM_score);

figure; hold on

for f = 1:numel(fcFields)
    fcName = fcFields{f};
    envFields = fieldnames(UP_DOWN_NSIM_score.(fcName));

    envRate = zeros(numel(envFields),1);
    diffScore = zeros(numel(envFields),1);

    for i = 1:numel(envFields)
        name = envFields{i};  
        envRate(i) = sscanf(name, 'env%d');
        diffScore(i) = UP_DOWN_NSIM_score.(fcName).(name);
    end

    [envRate, idx] = sort(envRate);
    diffScore = diffScore(idx);

    plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
         'DisplayName', fcName)
end

xlabel('Envelope rate')
ylabel('Differnece score')
title('Differnece score over envelope rate NSIM')
legend('Location','best')
grid on
hold off

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

    plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
         'DisplayName', fcName)
end

xlabel('Envelope rate')
ylabel('Differnece score')
title('Differnece score over envelope rate NSIM LOCAL')
legend('Location','best')
grid on
hold off

%% plot score over envRate NSIM_input
fcFields = fieldnames(UP_DOWN_NSIM_score_input);

figure; hold on

for f = 1:numel(fcFields)
    fcName = fcFields{f};
    envFields = fieldnames(UP_DOWN_NSIM_score_input.(fcName));

    envRate = zeros(numel(envFields),1);
    diffScore = zeros(numel(envFields),1);

    for i = 1:numel(envFields)
        name = envFields{i};  
        envRate(i) = sscanf(name, 'env%d');
        diffScore(i) = UP_DOWN_NSIM_score_input.(fcName).(name);
    end

    [envRate, idx] = sort(envRate);
    diffScore = diffScore(idx);

    plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
         'DisplayName', fcName)
end

xlabel('Envelope rate')
ylabel('Differnece score')
title('Differnece score over envelope rate NSIM input')
legend('Location','best')
grid on
hold off

% %% plot score over envRate NSIM_input_out_up
% fcFields = fieldnames(UP_DOWN_NSIM_score_input_out_up);
% 
% figure; hold on
% 
% for f = 1:numel(fcFields)
%     fcName = fcFields{f};
%     envFields = fieldnames(UP_DOWN_NSIM_score_input_out_up.(fcName));
% 
%     envRate = zeros(numel(envFields),1);
%     diffScore = zeros(numel(envFields),1);
% 
%     for i = 1:numel(envFields)
%         name = envFields{i};  
%         envRate(i) = sscanf(name, 'env%d');
%         diffScore(i) = UP_DOWN_NSIM_score_input_out_up.(fcName).(name);
%     end
% 
%     [envRate, idx] = sort(envRate);
%     diffScore = diffScore(idx);
% 
%     plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
%          'DisplayName', fcName)
% end
% 
% xlabel('Envelope rate')
% ylabel('Differnece score')
% title('Differnece score over envelope rate NSIM input out up')
% legend('Location','best')
% grid on
% hold off
% 
% %% plot score over envRate NSIM_input_out_down
% fcFields = fieldnames(UP_DOWN_NSIM_score_input_out_down);
% 
% figure; hold on
% 
% for f = 1:numel(fcFields)
%     fcName = fcFields{f};
%     envFields = fieldnames(UP_DOWN_NSIM_score_input_out_down.(fcName));
% 
%     envRate = zeros(numel(envFields),1);
%     diffScore = zeros(numel(envFields),1);
% 
%     for i = 1:numel(envFields)
%         name = envFields{i};  
%         envRate(i) = sscanf(name, 'env%d');
%         diffScore(i) = UP_DOWN_NSIM_score_input_out_down.(fcName).(name);
%     end
% 
%     [envRate, idx] = sort(envRate);
%     diffScore = diffScore(idx);
% 
%     plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
%          'DisplayName', fcName)
% end
% 
% xlabel('Envelope rate')
% ylabel('Differnece score')
% title('Differnece score over envelope rate NSIM input out down')
% legend('Location','best')
% grid on
% hold off

%% plot score over envRate SSIM
fcFields = fieldnames(UP_DOWN_SSIM_score);

figure; hold on

for f = 1:numel(fcFields)
    fcName = fcFields{f};
    envFields = fieldnames(UP_DOWN_SSIM_score.(fcName));

    envRate = zeros(numel(envFields),1);
    diffScore = zeros(numel(envFields),1);

    for i = 1:numel(envFields)
        name = envFields{i};  
        envRate(i) = sscanf(name, 'env%d');
        diffScore(i) = UP_DOWN_SSIM_score.(fcName).(name);
    end

    [envRate, idx] = sort(envRate);
    diffScore = diffScore(idx);

    plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
         'DisplayName', fcName)
end

xlabel('Envelope rate')
ylabel('Differnece score')
title('Differnece score over envelope rate SSIM')
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

    plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
         'DisplayName', fcName)
end

xlabel('Envelope rate')
ylabel('Differnece score')
title('Differnece score over envelope rate SSIM LOCAL')
legend('Location','best')
grid on
hold off

%% plot score over envRate SSIM_input
fcFields = fieldnames(UP_DOWN_SSIM_score_input);

figure; hold on

for f = 1:numel(fcFields)
    fcName = fcFields{f};
    envFields = fieldnames(UP_DOWN_SSIM_score_input.(fcName));

    envRate = zeros(numel(envFields),1);
    diffScore = zeros(numel(envFields),1);

    for i = 1:numel(envFields)
        name = envFields{i};  
        envRate(i) = sscanf(name, 'env%d');
        diffScore(i) = UP_DOWN_SSIM_score_input.(fcName).(name);
    end

    [envRate, idx] = sort(envRate);
    diffScore = diffScore(idx);

    plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
         'DisplayName', fcName)
end

xlabel('Envelope rate')
ylabel('Differnece score')
title('Differnece score over envelope rate SSIM input')
legend('Location','best')
grid on
hold off

% %% plot score over envRate SSIM_input_out_up
% fcFields = fieldnames(UP_DOWN_SSIM_score_input_out_up);
% 
% figure; hold on
% 
% for f = 1:numel(fcFields)
%     fcName = fcFields{f};
%     envFields = fieldnames(UP_DOWN_SSIM_score_input_out_up.(fcName));
% 
%     envRate = zeros(numel(envFields),1);
%     diffScore = zeros(numel(envFields),1);
% 
%     for i = 1:numel(envFields)
%         name = envFields{i};  
%         envRate(i) = sscanf(name, 'env%d');
%         diffScore(i) = UP_DOWN_SSIM_score_input_out_up.(fcName).(name);
%     end
% 
%     [envRate, idx] = sort(envRate);
%     diffScore = diffScore(idx);
% 
%     plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
%          'DisplayName', fcName)
% end
% 
% xlabel('Envelope rate')
% ylabel('Differnece score')
% title('Differnece score over envelope rate SSIM input out up')
% legend('Location','best')
% grid on
% hold off
% %% plot score over envRate SSIM_input
% fcFields = fieldnames(UP_DOWN_SSIM_score_input_out_down);
% 
% figure; hold on
% 
% for f = 1:numel(fcFields)
%     fcName = fcFields{f};
%     envFields = fieldnames(UP_DOWN_SSIM_score_input_out_down.(fcName));
% 
%     envRate = zeros(numel(envFields),1);
%     diffScore = zeros(numel(envFields),1);
% 
%     for i = 1:numel(envFields)
%         name = envFields{i};  
%         envRate(i) = sscanf(name, 'env%d');
%         diffScore(i) = UP_DOWN_SSIM_score_input_out_down.(fcName).(name);
%     end
% 
%     [envRate, idx] = sort(envRate);
%     diffScore = diffScore(idx);
% 
%     plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
%          'DisplayName', fcName)
% end
% 
% xlabel('Envelope rate')
% ylabel('Differnece score')
% title('Differnece score over envelope rate SSIM input out down')
% legend('Location','best')
% grid on
% hold off

%% plot score over envRate Eucl-Dist
fcFields = fieldnames(Channel_diff);

figure; hold on

for f = 1:numel(fcFields)
    fcName = fcFields{f};
    envFields = fieldnames(Channel_diff.(fcName));

    envRate = zeros(numel(envFields),1);
    diffScore = zeros(numel(envFields),1);

    for i = 1:numel(envFields)
        name = envFields{i};  
        envRate(i) = sscanf(name, 'env%d');
        diffScore(i) = Channel_diff.(fcName).(name);
    end

    [envRate, idx] = sort(envRate);
    diffScore = diffScore(idx);

    plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
         'DisplayName', fcName)
end

xlabel('Envelope rate')
ylabel('Difference score')
title('Difference score over envelope rate euclidic distance')
legend('Location','best')
grid on
hold off

%% plot score over envRate Eucl-Dist
fcFields = fieldnames(Channel_correlationDiff);

figure; hold on

for f = 1:numel(fcFields)
    fcName = fcFields{f};
    envFields = fieldnames(Channel_correlationDiff.(fcName));

    envRate = zeros(numel(envFields),1);
    diffScore = zeros(numel(envFields),1);

    for i = 1:numel(envFields)
        name = envFields{i};  
        envRate(i) = sscanf(name, 'env%d');
        diffScore(i) = Channel_correlationDiff.(fcName).(name);
    end

    [envRate, idx] = sort(envRate);
    diffScore = diffScore(idx);

    plot(envRate, diffScore, '-o', 'LineWidth', 2, ...
         'DisplayName', fcName)
end

xlabel('Envelope rate')
ylabel('Difference score')
title('Difference score over envelope rate Correlation')
legend('Location','best')
grid on
hold off

toc