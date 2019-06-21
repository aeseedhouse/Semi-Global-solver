function psi = evolve(t,tStep,M,L,psi0,calcH,threshold,varargin)
% Implements the method for solving the time dependent schrodinger equation
% detailed in arXiv:1611.06707
%
% psi = evolve(t,tStep,M,L,psi0,calcH,varargin)
% 
% Inputs:
% t - times at which to calculate psi(t)
% tStep - lengh of polynomial approximation interval for s_ext(t)
% M - number of terms to use in polynomial approximation of s_ext(t)
% L - number of terms to use in polynomial approximation of f_M(H0,t)
% psi0 - Dx1 vector containing initial state
% calcH - handle to a function which returns the hamiltonian
% varargin - any inputs recquired by calcH
% 
% Outputs:
% psi - Dxlength(t) matrix where psi(:,m) = psi(t(m))

    D = length(psi0);
    % timing
    nStep = ceil(t(end)/tStep); % number of approximation intervals
    x = -cos(pi*(0:M-1)/(M-1)); % chebyshev sampling points in domain [-1,1]
    tCheby = (x*tStep+tStep)/2; % sampling points within time interval
    tcalcPsi = [tCheby,tCheby(2:end)+tStep]; 
    % tCheby are times to calculate psi within a specific approximation interval 
    % tCheby(2:end)+tStep are times to calculate psi to use as the initial
    % guess for the next approximation interval

    psi = zeros(D,length(t));   % output - psi(:,m) = psi(t(m))
    psiGuess = psi0*ones(1,M);    % initial guess is just initial state

    % count = 0;  % count the number of times psi is recalculated

    for n=1:nStep
        tSample = tCheby+tStep*(n-1);
        % sample Hamiltonian
        H = calcH(tSample,varargin{:});
        H = -1i*H;  % taking hbar = 1;

        notConverged = 1;
    %     count = 0;  
        while notConverged
%             count = count + 1;
            [V,H0] = calcV(tCheby,psiGuess,H);
            psiNew = applyF(M,H0,tcalcPsi,V,min(D,L),threshold);
            if any(any(isnan(psiNew))) || any(any(isinf(psiNew)))
                error('diverged');
            end
            notConverged = norm(psiNew(:,M)-psiGuess(:,end))/norm(psiGuess(:,end)) > threshold;
            psiGuess = psiNew(:,1:M);
        end
        first = find(t >= tStep*(n-1),1,'first');   % index of first desired time stamp within current approximation interval
        last = find(t <= tStep*n,1,'last');   % index of last desired time stamp
        psi(:,first:last) = applyF(M,H0,t(first:last)-tStep*(n-1),V,min(D,L),threshold); % once converged, calculate solution for desired time stamps   
        psiGuess = psiNew(:,M:end);     % extrapolate current polynomial approximation into next time interval to provide initial guess

        % normalise uGuess in case of divergence
        for m=2:M
            psiGuess(:,m) = psiGuess(:,m)/norm(psiGuess(:,m));
        end
    end
end