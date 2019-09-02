function [gamma,vKrylov,w] = arnoldi(A,v,L)
% Arnoldi iteration to compute orthogonal basis for krylov space spanned by 
% v,A*v,A^2*v,...,A^(L-1)*v
% 
% Inputs:
% A - DxD matrix for which to calculate f(A)
% v - Dx1 state vector
% L - dimension of krylov space
% 
% Outputs:
% gamma - LxL matrix where gamma_i,j is the matrix element of A_i,j, in orthonormal basis of krylov space
% vMat - Dx(L+1) matrix where column n is nth orthonormal basis vector of krylov space
% w - representation of v in vMat basis

    gamma = zeros(L,L);
    D = size(v,1);
    vKrylov = zeros(D,L);
    
    % w = v in krylov basis
    w = zeros(D,1);
    w(1) = sqrt(v'*v); 
    if any(v) % if v is zero vector, return gamma, vKrylov and w as all zeros
        vKrylov(:,1) = v/w(1); % first basis vector is normalised input vector
        
        for jj = 1:L
            v_j1 = A*vKrylov(:,jj);
            
            % calculate projections of v_j1 onto krylov basis vectors and
            % remove from v_j1
            for ii=0:jj-1
                gamma(ii+1,jj) = vKrylov(:,ii+1)'*v_j1;
                v_j1 = v_j1 - gamma(ii+1,jj)*vKrylov(:,ii+1);
            end
            
            if jj < L   % don't save last vector and extra row of gamma
                if abs(v_j1) < 1e-6    % dimension of krylov space larger than dimension of v - gram schmidt starts producing zero vectors
                    gamma(jj+1,jj) = 0;
                    vKrylov(:,jj+1) = 0;
                else
                    gamma(jj+1,jj) = sqrt(v_j1'*v_j1);    % norm of v_j1
                    vKrylov(:,jj+1) = v_j1/gamma(jj+1,jj);    % normalise v_j1
                end
            end
        end
    end
end