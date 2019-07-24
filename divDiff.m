function cNewt = divDiff(fSample,xSample)
% Calculate expansion coefficients for newton interpolation of f(x) using
% samples of f at points x. Expansion coefficients are given by divided differences defined as:
% f[x1] = f(x1)
% f[x1,x2] = (f(x2)-f(x1))/(x2-x1)
% f[x1,...,xn] = (f[x2,...,xn]-f[x1,...,xn-1])/(xn-x1)
%
% f(x) ~= sum(m=0 to M-1) cNewt(:,m)*Rm(x)
% R0(x) = 1
% Rm(x) = product(n=1 to m) (x-xSample(n)) for m=1,2,...
% cNewt(:,m) = f[x1,..,xm]
% 
% Inputs:
% f - DxM matrix of sampled function values f(x) 
% x - DxM matrix of sampling points
% 
% Outputs:
% cNewt - DxM expansion coefficients for newton interpolation

    [D,M] = size(fSample);
    cNewt = zeros(D,M);
    diffs = fSample; % 0th order differences are sampled function values
    cNewt(:,1) = diffs(:,1);
    for d=2:M
        diffs = diff(diffs,1,2)./(xSample(:,d:M)-xSample(:,1:M-d+1)); % calculate (d-1)th order differences
        cNewt(:,d) = diffs(:,1);
    end
end