%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates and plots the outage probability for a 1D fluid antenna using
% both the average correlation model [16] and Jakes's model. 
% MAKES FIGURE 2 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------
% Initialization
%-------------------------------------------------------------------------
close all
clear
clc

addpath('Core/')

%-------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------
N = 100;                 % Number of ports 
W = 4;                   % Antenna size normalized by wavelength
U = 3;                   % Number of users

Nsamples = 1e5;          % Number of samples for Monte-Carlo simulation

gamdB = linspace(-5,15,50);     % SIR threshold (dB)
gam = 10.^(gamdB/10);           % SIR threshold (linear)

%-------------------------------------------------------------------------
% Jakes's model
%-------------------------------------------------------------------------
% Correlation matrix according to Jakes's model
Sigma_jakes = toeplitz(besselj(0,2*pi*(0:N-1)*W/(N-1)));

%-------------------------------------------------------------------------
% Average correlation model
%-------------------------------------------------------------------------
% [Eq. (5), 16]
mu = sqrt(2)*sqrt(hypergeom(0.5, [1 1.5], -pi^2*W^2) - ...
                        besselj(1, 2*pi*W)/(2*pi*W));

%---------------------------------------------------------------------
% Outage probabilities
%---------------------------------------------------------------------
% Jakes's outage probability
disp("Computing Jakes's outage probability...")
pout_jakes = SimOutage(Nsamples, gam, Sigma_jakes, U);

% Outage probability according to the average correlation model in [16]
disp("Computing average outage probability...")
pout_avg = CalcOutage(gam, N, mu^2, U, 'Quadrature',30);

%-------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------
figure(1)

semilogy(gamdB, pout_jakes, 'k--', 'linewidth', 2, 'DisplayName',...
                ['Jakes N=' num2str(N)]);

hold on; grid on;

semilogy(gamdB, pout_avg, 'k', 'linewidth', 2, 'DisplayName',...
                ['Average [16] N=' num2str(N)]);

l = legend;

set(gca, 'TickLabelInterpreter', 'latex','FontSize',18) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('$\gamma$ (dB)', 'FontSize', 18, 'Interpreter','latex');
ylabel('$P(SNR < \gamma)$', 'FontSize', 18, 'Interpreter','latex');
ylim([1e-3, 1])