function Q = calcqNewt(tNewt,dt)
% Calculate coefficients  to convert between expansion coefficients for newton interpolation and taylor interpolation
% s(t) ~= sum(m=0 to M-1) c_m*R_m(t) where R_m is the mth order newton basis polynomial
% and we want to convert to taylor-like expansion s(t) ~= sum(m=0 to M-1) s_m*t^m.

% s_m = sum(n=m to M-1) Q(n,m)*c_n

% Coefficients calculated assuming interpolation is over domain of length 4
% for numerical stability:   
% Q(0,0) = 1;
% Q(n+1,0) = -4*t(n)*q(n,0)/dt
% Q(n+1,m) = 4*(q(n,m-1)-t(n)*q(n,m)),  1<=m<=n
% Q(n+1,n+1) = 4*q(n,n)/dt      (eqs 220,226-228)

% Inputs:
% tNewt - newton sampling points in length 4 time interval
% dt - length of original time interval
% 
% Outputs:
% Q - MxM matrix of coefficients

    M = length(tNewt);
    factor = 4/dt;
    Q = zeros(M);
    Q(1, 1) = 1;    % First row:
    for n = 2:M
        Q(n, 1) = -factor*tNewt(n-1)*Q(n-1,1);
        for m = 2:n-1
            Q(n,m) = factor*(Q(n-1,m-1) - tNewt(n-1)*Q(n-1,m));
        end
        Q(n,n) = Q(n-1,n-1)*factor;
    end
end