function psi = applyF(M,H0,t,V,L,threshold)
% Calculate u = f_M(A,t)*v_M + sum(m=0 to M-1) v_m*t^m    (eqn 81)
% 
% f_M(H0,t) given by:
% f_M(H0,t) = (M!/H0^M)*(exp(H0*t)-sum(j=0 to M-1)(H0*t)^j/j!)  for A != 0
% f_M(H0,t) = t^M    for A = 0
% 
% Inputs:
% M - number of terms in polynomial approximation of s_ext(t)
% H0 - constant component of hamiltonian
% t - vector of times
% V - DxM matrix of where V(:,m) = v_m
% L - number of terms in polynomial approximation of f_M(H0,t)
% threshold - threshold for testing convergence in calcF
% 
% Outputs:
% vMat - LxL matrix where mth column is mth orthonormal basis vector for krylov space spanned by v,A*v,A^2*v,...,A^(L-1)*v
% gamma - LxL matrix where gamma_i,j is A_i,j in vMat basis 

    vM = V(:,end);
    timePowers = powers(t',M-1)';    % timePowers(ii,:) = t^(ii-1)
    
    if any(any(H0))
        [gamma,vKrylov,w] = arnoldi(H0,vM,L);
        eigs = eig(gamma);
        cap = capacity(mean(eigs),eigs);    % capacity of eigenvalue domain
        % calculate R_n(gamma)*w terms in newton approximation 
        % f(gamma,t) ~= sum(n=0 to M-1) a_n*R_n(gamma)*w
        Rnw = calcRn(w,gamma,eigs,L,cap); % newton basis polynomials       
        cNewt = divDiff(calcf(eigs,t,M,threshold).',eigs.'/cap,length(t),L).';
        f_vM = vKrylov*Rnw*cNewt;   % (eqn 189)
    else    % all zero hamiltonian
        f_vM = vM*t.^M;                
    end  
    psi = V(:,1:end-1)*timePowers + f_vM;
end

function cap = capacity(xp,x)
% calculate the capacity of the domain x defined by
% capacity = product(n=0 to N-1) |xp-xn|^(1/N)
    cap = 1;
    x = x(x~=xp); % remove any element of x equal to the test point xp;
    N = length(x);
    for n=1:N
        cap = cap*abs(x(n)-xp);
    end
    cap = cap^(1/N);
end

function f = calcf(z,t,M,threshold)
% adapted from Tal Ezer's implementation (supp materials arXiv:1611.06707)
% calculate f_M(z,t) (eqn. 82)
    f = ones(length(z),length(t));
    if size(t,2)==1
        t=t';
    end
    zt = z*t;
    direct = factorial(M)*eps./abs((z*t).^M) < threshold;
    tail = ~direct;
    f(direct) = exp(zt(direct));
    
    % Preliminary calculation for values which can be found using direct
    % application of eqn. 82
    for m=1:M
        f(direct) = m*(f(direct)-1)./zt(direct);
    end
    
    % Preliminary calculation of eqn. 89 for values which are too small to evaluate
    % using eqn. 82
    term = double(tail);
    j = 1;
    while max(max(tail))
        term(tail) = zt(tail).*term(tail)/(M+j);
        f(tail) = f(tail) + term(tail);
        j = j + 1;
        tail(tail) = abs(term(tail))./abs(f(tail)) > threshold;
    end
    f = f.*((ones(size(z))*t).^M);
end     

function Rn = calcRn(w,gamma,eigSamples,L,cap)
% Calculate the newton basis vectors Rn(z)*w in the approximation
% f(gamma)*w ~= sum(n=0 to L-1) an*Rn(z)*w
% R0(z) = 1;
% Rn(z) = product(m=1 to n) (z-zm)
% where zm's are sampling points (eigSamples)

    Rn = zeros(L,L);
    Rn(1, 1) = w(1);   % R0(gamma) = 1 => R0(gamma)*w = w;
    for n = 2:L
        Rn(:,n) = (gamma*Rn(:,n-1) - eigSamples(n-1)*Rn(:,n-1))/cap;
    end
end