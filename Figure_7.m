%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represents the OP achieved by slow-FAMA with a 1D fluid antenna
% under Jakes's correlation model and the two proposed approximations in
% Section III-B (Algorithm 1 and equal block sizes)
% MAKES FIG 7 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------
clc
close all
clear

addpath('Core/')

%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------
N = 120;                    % Number of ports

W = 8;                      % Antenna size (wavelength normalized)
U = 3;                      % Number of users

Nsamples = 1e6;             % Number of samples for Monte-Carlo simulation

gamdB = linspace(-5,15,21);    % SIR threshold (dB)
gam = 10.^(gamdB/10);          % SIR threshold (linear scale)

%-----------------------------------------------------------------
% Jakes's correlation
%----------------------------------------------------------------
Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N-1)*W/(N-1)));

rho = sort(eig(Sigma_jakes),'descend');

%----------------------------------------------------------------
% Block diagonal correlation matrix approximation
%-----------------------------------------------------------------
mu2 = 0.96;
Num_eig = sum(rho > N/100);

% Resulting block correlation when Algorithm 1 is applied
L = BlockCorrelation(N, rho, Num_eig, mu2);

% Resulting block correlation where all the blocks have equal size
Lequal = ones(Num_eig,1)*round(N/Num_eig);

%-----------------------------------------------------------------
% Outage Probabilities calculation/simulation
%----------------------------------------------------------------- 
% Jakes's outage probability
disp('Simulating Jakes outage probability...')
pout_jakes = SimOutage(Nsamples, gam, Sigma_jakes, U);

% Block outage probability -> Algorithm 1
disp('Calculating blocks outage probability...')
pout_blocks = CalcOutage(gam, L, mu2, U, 'Quadrature', 30);

% Block outage probability -> Equal Blocks
disp('Calculating blocks outage probability...')
pout_blocks_equal = CalcOutage(gam, Lequal, mu2, U, 'Quadrature', 30);

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

semilogy(gamdB, pout_jakes, 'k-o', 'linewidth', 2,'MarkerSize',8);
hold on; grid on;
semilogy(gamdB, pout_blocks, 'r-*', 'linewidth', 2,'MarkerSize',8);
semilogy(gamdB, pout_blocks_equal, 'b--d', 'linewidth', 2,'MarkerSize',8);

l = legend("Jakes'", "Algorithm 1", "Equal $L_b$");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('$\gamma$', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{out}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-5, 1])


