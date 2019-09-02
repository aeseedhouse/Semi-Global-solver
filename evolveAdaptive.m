function psi = evolveAdaptive(t,M,L,psi0,calcH,threshold,jumpThreshold,desiredStep,varargin)
% Implements the method for solving the time dependent schrodinger equation
% detailed in arXiv:1611.06707
%
% psi = evolveAdaptive(t,M,L,psi0,calcH,threshold,jumpThreshold,desiredStep,varargin)
% 
% Inputs:
% t - times at which to calculate psi(t)
% M - number of terms to use in polynomial approximation of s_ext(t)
% L - number of terms to used in applyF for polynomial approximation of f_M(H0,t) when
% calculating solution for psi(t)
% psi0 - Dx1 vector containing initial state
% calcH - handle to a function which calculates the hamiltonian at requested sampling times
% threshold - threshold for convergence checks
% jumpThreshold - threshold for maximum jump in Hamiltonian (leave empty if jump detection not needed)
% desiredStep - desired time step (leave empty if no desired step) 
% varargin - any inputs recquired by calcH
% 
% Outputs:
% psi - Dxlength(t) matrix where psi(:,m) = psi(t(m))


    % if no jumpThreshold specified, set it to infinity so it is never exceeded
    if isempty(jumpThreshold)
        jumpThreshold = inf;
    end
    
    D = length(psi0);
    
    % timing
    if ~isempty(desiredStep)
        tStep = 2*desiredStep;
    else
        tStep = t(end)/10;
    end
    x = -cos(pi*(0:M-1)/(M-1)); % chebyshev sampling points in domain [-1,1]
    tCheby = (x*tStep+tStep)/2; % sampling points within time step
    tCalc = [tCheby tCheby(2:end)+tStep]; % sampling points within current and next time step
    
    psi = zeros(D,length(t));   % psi(:,m) = psi(t(m))
    psiGuess = psi0*ones(1,M);    % guess for psi(t) in the first time step is the initial state
    
    % flags and indicators
    first = 1;  % index of first desired time stamp in current approximation interval     
    countMax = 20; % maximum number of iterations allowed when checking for convergence of solution    
    lastStep = 0; % last time step reached
    firstStep = 1; % first time step
    split = 0; % 1 if hamiltonian jump threshold exceeded or psi starting to diverge at current time step length
    prevSplit = 0; % 1 if previous time step was split, 0 otherwise
    changetStep = 0; % 1 if psi is diverging and time step needs to be reduced
%     tStep will be changed when:
%        - psi(t) at sampling points doesn't pass convergence checks
%        - elapsedTime + tStep runs over t(end), i.e. last time step has been reached
    
    elapsedTime = 0;
    Q = calcqNewt(tCheby,tCheby(end)-tCheby(1));  % newton interpolation conversion factors used in calcV (eqns 220,226-228)
    
    while elapsedTime < t(end)        
        if changetStep 
            tStep = newStep;
            tCheby = (x*tStep+tStep)/2; tCalc = [tCheby tCheby(2:end)+tStep]; % recalculate sampling times     
            
            % recalculate psiGuess for new sampling times
            if ~firstStep
%                 psiGuess(:,2:end) = applyF(M,H0,tCheby(2:end),V,min(D,L),threshold);
%                 psiGuess = psiGuess./vecnorm(psiGuess); 
                psiGuess = psiGuess(:,1)*ones(1,M);
            else
                psiGuess = psi0*ones(1,M); % if still in first time step, use psi0 as the guess
            end            
            Q = calcqNewt(tCheby,tCheby(end)-tCheby(1));    % recalculate newton interpolation conversion factors
            
            % reset flags
            changetStep = 0; 
%             split = 0; prevSplit = 0;
        end
        
        if elapsedTime + tStep > t(end) 
            % last time step has been reached
            lastStep = 1;
            tStep
            tStep = t(end)-elapsedTime;
            tCheby = (x*tStep+tStep)/2;     % recalculate sampling times
        end
        tSample = tCheby + elapsedTime;
        
        % sample Hamiltonian
        H = calcH(tSample,varargin{:}); % sample Hamiltonian for use in calcV
        H = -1i*H;  % H/(i*hbar) with hbar = 1;
        
        if abs(max(max(max(diff(H,1,3),[],3)))) > jumpThreshold && ~firstStep && ~prevSplit
            split = 1; % if maximum jump in hamiltonian exceeds threshold, split current time step
        else
            % calculate vj's
            notConverged = 1;
            currCount = 0;
            while notConverged
                currCount = currCount + 1;
                [V,H0] = calcV(tCheby,psiGuess,H,Q);
                psiNew = applyF(M,H0,tCalc,V,min(D,L),threshold);
                notConverged = norm(psiNew(:,M)-psiGuess(:,end))/norm(psiGuess(:,end)) > threshold;
                psiGuess = psiNew(:,1:M);
                if currCount > countMax
                    break;
                end
            end
        end
        
        if mean(abs(vecnorm(psiGuess)-1)) > threshold || split % divergence indicated by norm of psi increasing beyond 1+threshold
            change = 1;
            if (~firstStep || split)  && ~prevSplit 
                % non convergence may be due to large jump in hamiltonian,
                % try splitting up current time step
                jumpInd = findJump(H); % find when biggest jump in hamiltonian occurs
                % locate jump on finer grid                
                tStart = tSample(jumpInd); tEnd = tSample(jumpInd + 1);
                for ii=1:3
                    tFine = linspace(tStart,tEnd,100);
                    jumpInd = findJump(calcH(tFine,varargin{:}));
                    tStart = tFine(jumpInd); tEnd = tFine(jumpInd+1);
                end
                
                if jumpInd > 1 % if jumpInd = 1, there probably wasn't a jump - don't split the interval, reduce time step 
                    % first split interval goes from elapsedTime up to the
                    % time of the jump
                    if t(first) <= tFine(jumpInd) % if there is a desired time stamp in the first split interval
                        last = findtInd(first,elapsedTime,tFine(jumpInd)-elapsedTime,t);
                        tSplit = [t(first:last),tFine(jumpInd)];
                        psiSplit = evolveSplit(tSplit-elapsedTime,elapsedTime,M,L,psiGuess(:,1),calcH,threshold,(tFine(jumpInd)-elapsedTime)/4,varargin{:});
                        psi(:,first:last) = psiSplit(:,1:end-1);
                        first = last + 1;
                    else
                        psiSplit = evolveSplit(tFine(jumpInd)-elapsedTime,elapsedTime,M,L,psiGuess(:,1),calcH,threshold,(tFine(jumpInd)-elapsedTime)/4,varargin{:});
                    end                    
                    
                    % second split interval goes from just after time of jump to
                    % end of time step
                    if t(first) <= elapsedTime + tStep % if any desired time stamps are within the second split interval
                        last = findtInd(first,tFine(jumpInd),elapsedTime+tStep-tFine(jumpInd),t);                        
                        tSplit = [t(first:last),elapsedTime+tStep];
                        psiSplit = evolveSplit(tSplit-tFine(jumpInd+1),tFine(jumpInd+1),M,L,psiSplit(:,end),calcH,threshold,(elapsedTime+tStep-tFine(jumpInd+1))/4,varargin{:});
                        psi(:,first:last) = psiSplit(:,1:end-1);
                        first = last + 1;
                    else
                        psiSplit = evolveSplit(elapsedTime+tStep-tFine(jumpInd+1),tFine(jumpInd+1),M,L,psiSplit(:,end),calcH,threshold,(elapsedTime+tStep-tFine(jumpInd+1))/4,varargin{:});
                    end
                    
                    elapsedTime = elapsedTime + tStep;
                    psiGuess = psiSplit(:,end)*ones(1,M); % use psi(elapsedTime) as initial guess of psi(t) for next time step
                    prevSplit = 1; split = 0; change = 0; % set/reset flags
                end
            end
            
            if change % reduce time step if the time step could not be split
                changetStep = 1;    % reduce tStep
                if (~isempty(desiredStep)) && (mod(desiredStep,tStep) == 0)
                    % if psi starts diverging after tStep has been set to a
                    % factor of desiredStep, set tStep to the next lowest factor
                    newStep = desiredStep/(desiredStep/tStep+1);
                else
                    newStep = 0.99*tStep;
                end
            end
        elseif (~lastStep) && (~isempty(desiredStep)) && (mod(desiredStep,tStep) ~= 0)  % if solution starts converging but tStep is not a factor of desiredStep
            changetStep = 1;
            if tStep >= desiredStep
                newStep = desiredStep;
            else
                newStep = desiredStep/ceil(desiredStep/tStep);
            end            
        else % once converged, calculate solution for desired time stamps                        
            if t(first) <= elapsedTime+tStep  % if any desired time stamps are within current approximation interval                
                % find index of last desired time stamp
                last = findtInd(first,elapsedTime,tStep,t);
                psi(:,first:last) = applyF(M,H0,t(first:last)-elapsedTime,V,min(D,L),threshold);    % calculate psi at desired times            
                psi(:,first:last) = psi(:,first:last)./vecnorm(psi(:,first:last));  % normalise psi
                first = last + 1;
            end
            % extrapolate current polynomial approximation into next time interval to provide initial guess
            psiGuess = psiNew(:,M:end); 
            psiGuess = psiGuess./vecnorm(psiGuess); % normalise psiGuess in case of divergence
            
            elapsedTime = elapsedTime + tStep; 
            split = 0; prevSplit = 0; firstStep = 0; % reset flags
        end
    end
end



function jumpInd = findJump(H)
% find index of biggest jump in hamiltonian, H
    jumpInd = 1;
    maxJump = max(max(abs(H(:,:,2)-H(:,:,1))));
    for ii=2:size(H,3)-1
        jump = max(max(abs(H(:,:,ii+1)-H(:,:,ii))));
        if jump > maxJump
            jumpInd = ii;
            maxJump = jump;
        end
    end
end

function last = findtInd(first,elapsedTime,tStep,t)
% find index of last desired time to calculate psi(t) within current time step
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