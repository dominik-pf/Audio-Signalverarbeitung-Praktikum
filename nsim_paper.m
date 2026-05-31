function [nsim_val, nsim_map] = nsim_paper(A, B)

% A, B = neurograms / NAPs
% expected size: time x CF channel

A = double(A);
B = double(B);

if ~isequal(size(A), size(B))
    error("A and B must have the same size");
end

% Normalize to [0,1], like image similarity
A = A - min(A(:));
B = B - min(B(:));

A = A ./ max(A(:) + eps);
B = B ./ max(B(:) + eps);

% SSIM / NSIM constants
K1 = 0.01;
K2 = 0.03;
L  = 1;

C1 = (K1*L)^2;
C2 = (K2*L)^2;

% local window
win = fspecial("gaussian", [11 11], 1.5);

% local means
muA = imfilter(A, win, "replicate");
muB = imfilter(B, win, "replicate");

% local variances and covariance
sigmaA2 = imfilter(A.^2, win, "replicate") - muA.^2;
sigmaB2 = imfilter(B.^2, win, "replicate") - muB.^2;
sigmaAB = imfilter(A.*B, win, "replicate") - muA.*muB;

% NSIM map
nsim_map = ((2*muA.*muB + C1) .* (2*sigmaAB + C2)) ./ ((muA.^2 + muB.^2 + C1) .* (sigmaA2 + sigmaB2 + C2));

% scalar NSIM score
nsim_val = mean(nsim_map(:), "omitnan");

end