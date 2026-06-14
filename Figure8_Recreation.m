%% Continuous heatmap: x = Fc, y = envRate
% color = NSIM local difference
% black = not tested because envRate > Fc

fcFields = fieldnames(UP_DOWN_NSIM_score_local);

%% numeric Fc values
FcVals = zeros(numel(fcFields),1);

for f = 1:numel(fcFields)
    FcVals(f) = sscanf(fcFields{f}, 'Fc%d');
end

[FcVals, sortIdx] = sort(FcVals);
fcFields = fcFields(sortIdx);

%% envRate axis: your generation only goes up to 1500 Hz
maxEnv = 1500;
envGrid = linspace(0, maxEnv, 500);

M = NaN(numel(envGrid), numel(FcVals));
% rows = envRate
% cols = Fc

for f = 1:numel(fcFields)

    fcName = fcFields{f};
    Fc = FcVals(f);

    envFields = fieldnames(UP_DOWN_NSIM_score_local.(fcName));

    envRates = zeros(numel(envFields),1);
    scores   = zeros(numel(envFields),1);

    for e = 1:numel(envFields)
        envName = envFields{e};

        envRates(e) = sscanf(envName, 'env%d');
        scores(e)   = UP_DOWN_NSIM_score_local.(fcName).(envName);
    end

    [envRates, idx] = sort(envRates);
    scores = scores(idx);

    %% interpolation only in valid region
    valid = envGrid >= min(envRates) & ...
            envGrid <= max(envRates) & ...
            envGrid <= Fc;

    M(valid,f) = interp1(envRates, scores, envGrid(valid), 'linear');

    %% optional: extend constant from last measured envRate up to limit
    fillTo = min(Fc, maxEnv);

    extend = envGrid > max(envRates) & ...
             envGrid <= fillTo;

    M(extend,f) = scores(end);

    %% impossible region stays black
    M(envGrid > Fc, f) = NaN;
end

threshold = [0.1 0.2 0.4 0.5 0.6 0.7];
crossEnv = NaN(length(threshold), numel(FcVals));
for k1 = 1:length(threshold)
    for z=1:numel(FcVals)
        col = M(:,z);
        indx05 = find(col>threshold(k1), 1, 'first');
        if ~isempty(indx05)
            crossEnv(k1,z) = envGrid(indx05);
        end
    end
end

%% Plot
figure;

h = imagesc(FcVals, envGrid, M);
axis xy;

set(gca, 'Color', 'k');
set(h, 'AlphaData', ~isnan(M));

xlabel('Fc (Hz)');
ylabel('Envelope rate (Hz)');
title('NSIM local difference score over Fc and envelope rate');

cb = colorbar;
ylabel(cb, 'Difference score');

colormap(parula);

validVals = M(~isnan(M));
caxis([min(validVals) max(validVals)]);

set(gca, 'XScale', 'log');
xticks(FcVals);
xticklabels(string(FcVals));

ylim([0 maxEnv]);

hold on
for k=1:length(threshold)
    % threshold points
    plot(FcVals, crossEnv(k,:), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'LineWidth', 1.5);
    
    % line through valid points
    validCross = ~isnan(crossEnv(k,:));
    
    plot(FcVals(validCross), crossEnv(k,validCross), 'r-', 'LineWidth', 2);

    % Label at the end of the line
    last = find(validCross, 1, 'last');
    if ~isempty(last)
        text(FcVals(last), crossEnv(k,last)+20, sprintf('%.1f', threshold(k)),'color','r', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end
hold off
grid on