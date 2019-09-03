function H = rabiHamOffRes(t,varargin)
% rotating frame hamiltonian for ESR pulse to rotate single spin about X
    fRabi = varargin{1};
    
    H = zeros(2,2,length(t));
    for ii=1:length(t)
        H(1,2,ii) = 2*pi*fRabi/2 + 2*pi*fRabi/2*cos(t(1,ii)*2*pi*fRabi);    % taking hbar = 1
        H(2,1,ii) = 2*pi*fRabi/2 + 2*pi*fRabi/2*cos(t(1,ii)*2*pi*fRabi);
    end
end
    
