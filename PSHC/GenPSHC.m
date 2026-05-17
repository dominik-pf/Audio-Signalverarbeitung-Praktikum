function [out, t] = GenPSHC(f0, M, N, k, dur, fs, direction)
    %function [out, t] = GenPSHC(f0, M, N, k, dur, fs);
    % generates a pulse-spreading harmonic complex (PSHC) following Hilkhuysen,
    % Gaston & Macherey, Olivier (201X) "Optimizing Pulse-Spreading Harmonic 
    % Complexes to minimize intrinsic modulations after auditory filtering". JASA.
    %
    % where:
    %  out  = vector representing PSHC in time domain
    %  t    = time vector corresponding to out                 [s]
    %
    %  f0   = fundamental frequency                           [Hz]
    %  M    = lowest harmonic # to include in PSHC
    %  N    = highest harmonic # to include in PSHC       
    %  k    = order of PSHC
    %  dur  = duration of signal                               [s]
    %  fs   = sampling rate                                   [Hz]
    
    t         = [1:dur*fs]/fs;         % time vector
    out       = zeros(1, length(t));   % initialize output vector
    HarmRank  = [M:N];                 % harmonics in PSHC
    j         = mod([M:N], k) + 1;     % eq(3)
    % r         = randperm(k)           % randomize the order of the different subcomplexes
    if direction == "down"
        r = j;
    elseif direction == "up"
        r = k - j+1;
    end
    u         = unifrnd(0,1,k,1);      % randomize the relative phase relation between components within each subcomplex
   
    for index = 1:length(HarmRank)                           % eq(1) summation
        % original but diviates from paper
        % phi= 2*pi()*((HarmRank(index)/k^2)*r(j(index)) +  ... % eq(2) phase relation
        %               u(j(index)));
        phi= 2*pi()*((HarmRank(index)/k^2)*r(index) ); % eq(2) phase relation

        out = out                                      +  ...
              sin(2*pi()*f0*HarmRank(index)*t  + phi);        % eq(1) harmonic series
    end                  

end

% By Pete D. Best
% 04/29/14