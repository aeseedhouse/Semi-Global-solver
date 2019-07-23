function [V,H0] = calcV(tSample,uGuess,H,Q) 
% Calculate vj's in:
% u(t) ~= fM(G0,t)v_M + sum(j=0 to M-1) t^j*v_j   (eqn 81)
% 
% Inputs:
% tSample - length D vector of sampling times
% uGuess - DxM matrix where column m is the initial guess for u(t) at time tSample(m)
% H - DxDxM matrix containing hamiltonian at times tSample
% 
% Outputs:
% V - Dx(M+1) matrix where column m is v_m
% H0 - constant part of hamiltonian
% H(t) = H0 + H'(t)

    D = size(uGuess,1); % dimension of qubit system
    M = length(tSample);

    tNewt = tSample*4/(tSample(end)-tSample(1)); % convert sampling times to length 4 domain for stability of newton interpolation
    H0 = H(:,:,round(M/2)); % constant part of hamiltonian

    %% calculate s_ext(t) at sampling points
    % s_ext(t) = H'(t)*u(t)
    %          = H(t)*u(t) - H0*u(t)  (eqn 58,59)
    s_ext = zeros(D,M);
    for m=1:M
        s_ext(:,m) = H(:,:,m)*uGuess(:,m);
    end
    s_ext = s_ext - H0*uGuess;
    
    %% Newton interpolation of s_ext(t)
    % s_ext(t) ~= sum(m=0 to M-1) cm*Rm(t)
    % Rm(t) - newton basis polynomial, cm = cNewt(:,m)
    cNewt = divDiff(s_ext,tNewt,D,M);
    
    % convert cm's to sm's  (eqn 226-228)
%     Q = calcqNewt(tSample,tSample(end)-tSample(1)); 
    sTaylor = cNewt*Q;

    % convert sm's to vm's    (eqn 83)
    V = zeros(D,M+1);
    V(:,1) = uGuess(:,1);
    for m=2:M+1
        V(:,m) = (H0*V(:,m-1) + sTaylor(:,m-1))/(m-1);
    end
end