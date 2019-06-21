clear;clc;close all;

%% initialise
% parameters
psi0 = [1;0]; % initial state
D = length(psi0);
M = 7;  % number of points to sample in polynomial series approximation of s_ext
L = M;    % number of points to sample in polynomial approximation f(G,t)
fRabi = 1e6;    % rabi frequency
threshold = 1e-6;   % threshold for convergence checks

t = linspace(0,1,100)*1.5e-6;     % desired times to calculate psi(t)
tStep = 1e-6;   % length of domain for approximation of s_ext(t)
calcH = @rabiHam;   % handle for function which returns the hamiltonian at given times

%% calculate psi(t)
psi = evolve(t,tStep,M,L,psi0,calcH,threshold,fRabi);

%% Plot
pX = zeros(size(t));
pY = pX;
pZ = pX;
for ii=1:length(t)
    pZ(ii) = abs(psi(1,ii))^2;
    pX(ii) = abs([1/sqrt(2),1/sqrt(2)]*psi(:,ii))^2;
    pY(ii) = abs([1/sqrt(2),-1i/sqrt(2)]*psi(:,ii))^2;
end

plot(t,pX)
hold on
plot(t,pY)
plot(t,pZ)
xlabel('Time')
ylabel('Probability')
legend('p(X)','p(Y)','p(Z)')
ylim([0 1])