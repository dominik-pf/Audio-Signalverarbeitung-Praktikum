%% Create PSHC
function x = generate_pshc(t, f0, harmonics, k, direction)

    x = zeros(size(t));


    for i = harmonics

        j = mod(i, k) + 1;

        if direction == "down"
            rj = j;
        elseif direction == "up"
            rj = k - j + 1;
        else
            error("direction must be 'down' or 'up'")
        end

        phi = 2*pi*rj*i/(k^2);

        x = x + sin(2*pi*f0*i*t + phi);
    end
end