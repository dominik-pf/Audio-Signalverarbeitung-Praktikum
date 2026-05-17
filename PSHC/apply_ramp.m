%% Ramp function
function y = apply_ramp(x, fs, rampDur)

    nRamp = round(rampDur * fs);
    ramp = 0.5 - 0.5*cos(pi*(0:nRamp-1)/nRamp);

    env = ones(size(x));
    env(1:nRamp) = ramp;
    env(end-nRamp+1:end) = fliplr(ramp);

    y = x .* env;
end
