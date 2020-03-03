function [psi,H0,V] = evolveSplit(t,tOffset,M,L,psi0,calcH,threshold,varargin)
% Same as evolveAdaptive but without jump detection or splitting time steps
% detailed in arXiv:1611.06707
%
% psi = evolve(t,tStep,M,L,psi0,calcH,varargin)
% 
% Inputs:
% t - times at which to calculate psi(t)
% M - number of terms to use in polynomial approximation of s_ext(t)
% L - number of terms to used in applyF for polynomial approximation of f_M(H0,t) when
% calculating solution for psi(t)
% psi0 - Dx1 vector containing initial state
% calcH - handle to a function which calculates the hamiltonian at requested
% sampling times
% threshold - threshold for convergence checks
% jumpThreshold - threshold for maximum jump in Hamiltonian (leave empty if jump detection not needed)
% desiredStep - desired time step (leave empty if no desired step) 
% varargin - any inputs recquired by calcH
% 
% Outputs:
% psi - Dxlength(t) matrix where psi(:,m) = psi(t(m))

    D = length(psi0);
    Mfactorial = factorial(M);
    % timing
    desiredStep = [];
    tStep = t(end);
    x = -cos(pi*(0:M-1)/(M-1)); % chebyshev sampling points in domain [-1,1]
    tCheby = (x*tStep+tStep)/2; % sampling points within time interval
    tCalc = [tCheby tCheby(2:end)+tStep];
    
    psi = zeros(D,length(t));   % psi(:,m) = psi(t(m))
    psiGuess = psi0*ones(1,M);    % guess for psi(t) in the first time step is the initial state
    
    % flags and indicators
    first = 1;  % index of first desired time stamp in current approximation interval  
    countMax = 20; % maximum number of iterations allowed when checking for convergence of solution
    changetStep = 0;
    lastStep = 0;
%     tStep will be changed when:
%        - psi(t) at sampling points doesn't pass convergence checks
%        - elapsedTime + tStep runs over t(end), i.e. last time step has been reached
    
    elapsedTime = 0;
    Q = calcqNewt(tCheby,tCheby(end)-tCheby(1));  % newton interpolation conversion factors used in calcV (eqns 220,226-228)
    
    psi(:,t==0) = psi0.*ones(D,sum(t==0)); 
    
    while elapsedTime < t(end)        
        if changetStep 
            tStep = newStep;
            tCheby = (x*tStep+tStep)/2;     % recalculate sampling times 
            tCalc = [tCheby tCheby(2:end)+tStep];
            % recalculate psiGuess for new sampling times
            if elapsedTime > 0
                psiGuess(:,2:end) = applyF(M,Mfactorial,H0,tCheby(2:end),V,min(D,L),threshold);
                psiGuess = psiGuess./vecnorm(psiGuess);
            else
                psiGuess = psi0*ones(1,M);
            end
            
            Q = calcqNewt(tCheby,tCheby(end)-tCheby(1));    % recalculate newton interpolation conversion factors
            
            % reset flags
            changetStep = 0;
        end
        
        if elapsedTime + tStep > t(end) 
            % last time step has been reached
            lastStep = 1;
            tStep = t(end)-elapsedTime;
            tCheby = (x*tStep+tStep)/2;     % recalculate sampling times
        end
        tSample = tCheby + elapsedTime;
        
        % sample Hamiltonian
        H = calcH(tSample+tOffset,varargin{:});
        H = -1i*H;  % H/(i*hbar) with hbar = 1;

        % calculate vj's 
        notConverged = 1;
        currCount = 0;  
        while notConverged
            currCount = currCount + 1;
            [V,H0] = calcV(tCheby,psiGuess,H,Q);
            psiNew = applyF(M,Mfactorial,H0,tCalc,V,min(D,L),threshold);
            notConverged = norm(psiNew(:,M)-psiGuess(:,end))/norm(psiGuess(:,end)) > threshold;
            psiGuess = psiNew(:,1:M);
            if currCount > countMax
                break;
            end
        end
        if mean(abs(vecnorm(psiGuess)-1)) > threshold % divergence indicated by norm of psi increasing beyond 1+threshold
            changetStep = 1;    % reduce tStep
            if (~isempty(desiredStep)) && (mod(desiredStep,tStep) == 0)
                % if psi starts diverging after tStep has been set to a
                % factor of desiredStep, set tStep to the next lowest factor
                newStep = desiredStep/(desiredStep/tStep+1);
            else
                newStep = 0.99*tStep;
            end
        elseif (~lastStep) && (~isempty(desiredStep)) && (mod(desiredStep,tStep) ~= 0)  % if solution starts to converge but tStep is not a factor of desiredStep
            changetStep = 1;
            if tStep >= desiredStep
                newStep = desiredStep;
            else
                newStep = desiredStep/ceil(desiredStep/tStep);
            end            
        else
            % once converged, calculate solution for desired time stamps
            if t(first) <= elapsedTime+tStep  % if any desired time stamps are within current approximation interval
                
                % find index of last desired time stamp
                last = findtInd(first,elapsedTime,tStep,t);
                psi(:,first:last) = applyF(M,Mfactorial,H0,t(first:last)-elapsedTime,V,min(D,L),threshold);    % calculate psi at desired times            
                psi(:,first:last) = psi(:,first:last)./vecnorm(psi(:,first:last));  % normalise psi
                
                first = last + 1;
            end
            % extrapolate current polynomial approximation into next time interval to provide initial guess
            psiGuess = psiNew(:,M:end);
            
            % normalise psiGuess in case of divergence
            psiGuess = psiGuess./vecnorm(psiGuess);
            elapsedTime = elapsedTime + tStep;
        end
    end
end

function last = findtInd(first,elapsedTime,tStep,t)
    last = first + 1;
    if elapsedTime + tStep == t(end)
        last = length(t);
    else
        while t(last) <= elapsedTime + tStep
            last = last + 1;
        end
        last = last - 1;
    end
end