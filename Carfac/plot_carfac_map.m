function h = plot_carfac_map(t, CFs, data, hin)

if nargin > 3
    h = hin;
else
    h = axes;
end

axes(h)

if size(data,1) ~= length(t)
    data = data.';
end

imagesc(t, log10(CFs/1e3), data.')
axis xy

yticks = [0.125 0.5 2 8 16];
set(gca,'YTick',log10(yticks))
set(gca,'YTickLabel',yticks)

xlabel('Time (s)')
ylabel('CF (kHz)')
hcb = colorbar;
set(get(hcb,'ylabel'),'string','Amplitude')
end