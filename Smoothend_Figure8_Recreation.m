%% Continuous heatmap from smoothed fitted functions

fcFields = fieldnames(UP_DOWN_NSIM_fitting);

%% numeric Fc values
FcVals = zeros(numel(fcFields),1);

for f = 1:numel(fcFields)
    FcVals(f) = sscanf(fcFields{f}, 'Fc%d');
end

[FcVals, sortIdx] = sort(FcVals);
fcFields = fcFields(sortIdx);

%% log envelope grid
envMin = 1;
maxEnv = 1500;
envGrid = logspace(log10(envMin), log10(maxEnv), 500)';

M = NaN(numel(envGrid), numel(FcVals));

for f = 1:numel(fcFields)

    fcName = fcFields{f};
    Fc = FcVals(f);

    envFine = UP_DOWN_NSIM_fitting.(fcName).envFine(:);
    diffFit = UP_DOWN_NSIM_fitting.(fcName).diffFit(:);

    valid = envGrid >= min(envFine) & ...
            envGrid <= max(envFine) & ...
            envGrid <= Fc;

    M(valid,f) = interp1(envFine, diffFit, envGrid(valid), 'pchip');

    extend = envGrid > max(envFine) & ...
             envGrid <= min(Fc,maxEnv);

    M(extend,f) = diffFit(end);

    M(envGrid > Fc,f) = NaN;
end

%% threshold crossings
threshold = [0.1 0.2 0.4 0.5 0.6 0.7];
crossEnv = NaN(length(threshold), numel(FcVals));

for k = 1:length(threshold)
    for z = 1:numel(FcVals)
        col = M(:,z);
        idx = find(col > threshold(k), 1, 'first');

        if ~isempty(idx)
            crossEnv(k,z) = envGrid(idx);
        end
    end
end

%% pixel edges for centered log pixels
xCenters = FcVals(:)';
xEdges = zeros(1,numel(xCenters)+1);
xEdges(2:end-1) = sqrt(xCenters(1:end-1).*xCenters(2:end));
xEdges(1) = xCenters(1)^2 / xEdges(2);
xEdges(end) = xCenters(end)^2 / xEdges(end-1);

yCenters = envGrid(:);
yEdges = zeros(numel(yCenters)+1,1);
yEdges(2:end-1) = sqrt(yCenters(1:end-1).*yCenters(2:end));
yEdges(1) = yCenters(1)^2 / yEdges(2);
yEdges(end) = yCenters(end)^2 / yEdges(end-1);

%% Plot
figure;
                                                                                xCenters = 1:length(FcVals);
                                                                                xEdges = 0.5:(length(FcVals)+0.5);

                                                                                [Xe,Ye] = meshgrid(xEdges,yEdges);

                                                                                C = NaN(numel(yEdges),numel(xEdges));
                                                                                C(1:end-1,1:end-1) = M;
                                                                                
                                                                                h = surf(Xe,Ye,zeros(size(C)),C);
                                                                                view(2)
                                                                                set(h,'EdgeColor','none')
                                                                                grid off
                                                                                set(gca,'XScale','linear')


% [Xe, Ye] = meshgrid(xEdges, yEdges);
% 
% C = NaN(numel(yEdges), numel(xEdges));
% C(1:end-1, 1:end-1) = M;
% 
% h = surf(Xe, Ye, zeros(size(C)), C);
% view(2);
% set(h, 'EdgeColor', 'none');

axis xy;
set(gca, 'Color', 'k');
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');


xlabel('Fc (Hz)');
ylabel('Envelope rate (Hz)');
title('NSIM local difference score over Fc and envelope rate');

cb = colorbar;
ylabel(cb, 'Difference score');

colormap(parula);

validVals = M(~isnan(M));
caxis([min(validVals) max(validVals)]);

xticks(FcVals);
xticklabels(string(FcVals));

ylim([envMin maxEnv]);

% hold on
% for k = 1:length(threshold)
% 
%     validCross = ~isnan(crossEnv(k,:));
% 
%     plot(FcVals(validCross), crossEnv(k,validCross), ...
%          'r-', 'LineWidth', 2);
% 
%     plot(FcVals(validCross), crossEnv(k,validCross), ...
%          'ro', 'MarkerFaceColor', 'r', ...
%          'MarkerSize', 7, 'LineWidth', 1.5);
% 
%     last = find(validCross, 1, 'last');
%     if ~isempty(last)
%         text(FcVals(last)*1.05, crossEnv(k,last)*1.08, ...
%              sprintf('%.1f', threshold(k)), ...
%              'Color','r', 'FontSize',10, ...
%              'HorizontalAlignment','left', ...
%              'VerticalAlignment','middle');
%     end
% end
% hold off
% grid on

hold on

xPlot = 1:length(FcVals);

for k = 1:length(threshold)

    validCross = ~isnan(crossEnv(k,:));

    plot(xPlot(validCross), crossEnv(k,validCross), ...
         'r-', 'LineWidth', 2);

    plot(xPlot(validCross), crossEnv(k,validCross), ...
         'ro', 'MarkerFaceColor', 'r', ...
         'MarkerSize', 7, 'LineWidth', 1.5);

    last = find(validCross, 1, 'last');
    if ~isempty(last)
        text(xPlot(last)+0.15, crossEnv(k,last)*1.08, ...
             sprintf('%.1f', threshold(k)), ...
             'Color','r', 'FontSize',10, ...
             'HorizontalAlignment','left', ...
             'VerticalAlignment','middle');
    end
end

hold off
grid off

xticks(1:length(FcVals))
xticklabels(string(FcVals))
xlabel('Fc (Hz)')