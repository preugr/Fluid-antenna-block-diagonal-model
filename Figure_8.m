%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the evolution of OP as a 1D fluid antenna is saturated (i.e., more
% ports are added while keeping the aperture constant). Jakes's correlation
% model is compared with the proposed block-diagonal approximation and the
% average correlation model in [16]
% MAKES FIGURE 8 IN THE PAPER
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
Nd = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30]; % Port density per wavelength.

W = 6;                       % Antenna size (wavelength normalized)
U = [3 4 5];                 % Number of users

Nsamples = 2e5;              % Number of samples for Monte-Carlo simulation

gam = 1;                     % SIR threshold (linear scale)

%-------------------------------------------------------------------------
% Pre-allocation
%-------------------------------------------------------------------------
% Jakes's model
pout_jakes = zeros(length(U),length(Nd));

% Block-diagonal correlation model 
pout_blocks = zeros(length(U),length(Nd));

% Average correlation model in [16]
pout_avg = zeros(length(U),length(Nd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP OVER FLUID ANTENNA NUMBER OF PORTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ksize = 1:length(Nd)

    % User feedback
    disp(['Iter ' num2str(ksize) ' out of ' num2str(length(Nd))]);

    %-----------------------------------------------------------------
    % Jake's correlation
    %----------------------------------------------------------------
    N = W * Nd(ksize);
    Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N-1)*W/(N-1)));

    rho = sort(eig(Sigma_jakes),'descend');

    %----------------------------------------------------------------
    % Block diagonal correlation matrix approximation -> Algorithm 1
    %-----------------------------------------------------------------
    mu2 = 0.97;
    Num_eig = sum(rho > N/100);
    L = BlockCorrelation(N, rho, Num_eig, mu2);

    %----------------------------------------------------------------
    % Average correlation factor -> [16]
    %-----------------------------------------------------------------
    mu_avg = sqrt(2)*sqrt(hypergeom(0.5, [1 1.5], -pi^2*W^2) - ...
                    besselj(1, 2*pi*W)/(2*pi*W));

    % Loop over number of users
    for ku = 1:length(U)
    
        %-----------------------------------------------------------------
        % Outage Probabilities calculation/simulation
        %----------------------------------------------------------------- 
        % Jakes's correlation
        pout_jakes(ku,ksize) = SimOutage(Nsamples, gam, Sigma_jakes, U(ku));
    
        % Block outage probability
        pout_blocks(ku,ksize) = CalcOutage(gam, L, mu2, U(ku), 'Quadrature', 25);
    
        % Average correlation model  
        pout_avg(ku,ksize) = CalcOutage(gam, N, mu_avg^2, U(ku),'Quadrature',25);
    
    end
end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for ku = length(U):-1:1

    semilogy(Nd, pout_jakes(ku,:), 'k-*', 'linewidth', 2);
    hold on; grid on;
    semilogy(Nd, pout_blocks(ku,:), 'b-o', 'linewidth', 2);
    semilogy(Nd, pout_avg(ku,:), '--v', 'Color', [0.4660 0.6740 0.1880], 'linewidth', 2);
end

l = legend("Jake's", "Blocks (29)", "Constant [16]");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 20, 'Interpreter','latex');
xlabel('$N/W$', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{out}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-4, 1])
xlim([Nd(1) Nd(end)])



