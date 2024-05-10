%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and plots the OP in terms of the SIR threshold for a 1D fluid
% antenna under different correlation models: 1) Jakes's, 2) Proposed block
% diagonal approximation (Eq. (29)), 3) Approximated block-diagonal 
% (Eq. (35)), 4) Average constant model in [16], 5) i.i.d. (Eq. (38))
% MAKES FIGURE 10 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
N = 100;                     % Number of ports

W = 5;                       % Antenna size (wavelength normalized)
U = [3 5];                   % Number of users

Nsamples = 5e5;              % Number of samples for Monte-Carlo simulation

gamdB = linspace(-10,10,20);     % SIR threshold (dB)
gam = 10.^(gamdB/10);            % SIR threshold (linear scale)

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------
% Jakes's model
pout_jakes = zeros(length(U),length(gam));

% Block-diagonal approximation (Lemma 2)
pout_blocks = zeros(length(U),length(gam));

% Approximated block-diagonal (Corollary 2)
pout_blocks_approx = zeros(length(U),length(gam));

% i.i.d. system
pout_iid = zeros(length(U),length(gam));

% Average correlation model in [16]
pout_avg = zeros(length(U),length(gam));

%-----------------------------------------------------------------
% Jake's correlation
%----------------------------------------------------------------
Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N-1)*W/(N-1)));

rho = sort(eig(Sigma_jakes),'descend');

%----------------------------------------------------------------
% Block diagonal correlation matrix approximation
%-----------------------------------------------------------------
mu2 = 0.97;
Num_eig = sum(rho > N/100);

% Algorithm 1
L = BlockCorrelation(N, rho, Num_eig, mu2);

% Average correlation coefficient in [16]
mu_avg = sqrt(2)*sqrt(hypergeom(0.5, [1 1.5], -pi^2*W^2) - ...
                besselj(1, 2*pi*W)/(2*pi*W));

% Loop over number of users
for ku = 1:length(U)

    % User feedback
    disp(['Iter ' num2str(ku) ' out of ' num2str(length(U))]);

    %-----------------------------------------------------------------
    % Outage Probabilities calculation/simulation
    %----------------------------------------------------------------- 
    % Jakes's outage probability
    pout_jakes(ku,:) = SimOutage(Nsamples, gam, Sigma_jakes, U(ku));

    % Block outage probability -> Evaluated through quadrature expression
    % in Eq. (31)
    pout_blocks(ku,:) = CalcOutage(gam, L, mu2, U(ku), 'Quadrature', 30);

    % Average correlation model in [16]  
    pout_avg(ku,:) = CalcOutage(gam, N, mu_avg^2, U(ku), 'Quadrature',30);

    % I.I.D. outage
    pout_iid(ku,:) = (1 - 1./(1+gam).^(U(ku)-1)).^Num_eig;

    % Approximated blocks outage probability -> Evaluated through
    % quadrature expression in Eq. (37)
    pout_blocks_approx(ku,:) = ApproxOutage(gam, mu2, U(ku), L, 'Quadrature', 30);
end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for ku = length(U):-1:1

    semilogy(gamdB, pout_jakes(ku,:), 'k-*', 'linewidth', 2);
    hold on; grid on;
    semilogy(gamdB, pout_blocks(ku,:), 'b-o', 'linewidth', 2);
    semilogy(gamdB, pout_avg(ku,:), '--v', 'Color', [0.4660 0.6740 0.1880], 'linewidth', 2);
    semilogy(gamdB, pout_iid(ku,:), 'r.-','Marker','square', 'linewidth', 2);
    semilogy(gamdB, pout_blocks_approx(ku,:), 'b--^', 'linewidth', 2);
end

l = legend("Jake's", "Blocks Eq. (29)", "Constant [12]", "i.i.d. Eq. (37)", "Blocks Eq. (35)");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('$\gamma$ (dB)', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{out}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-5, 1])



