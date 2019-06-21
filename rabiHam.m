function H = rabiHam(t,varargin)
% rotating frame hamiltonian for ESR pulse to rotate single spin about X
    fRabi = varargin{1};
    
    H = zeros(2,2,length(t));
    H(1,2,:) = 2*pi*fRabi/2;    % taking hbar = 1
    H(2,1,:) = 2*pi*fRabi/2;
end
    
