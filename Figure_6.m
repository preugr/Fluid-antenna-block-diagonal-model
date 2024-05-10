%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates and represent the gain, in terms of OP, achieved by each block
% of the proposed correlation model, i.e., Eq. (39).
% MAKES FIG 6 IN THE PAPER
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
gam = [0.5 1];                          % SIR threshold
mu2 = 0.97;                             % mu^2
U = [4 5 6 7];                          % Number of users 
L = 2:100;                              % Block size

%-------------------------------------------------------------------------
% Delta Outage Probability
%-------------------------------------------------------------------------
% Pre-allocating
DeltaPout = zeros(length(U), length(L), length(gam));

% Loop over SIR thresholds
for kg = 1:length(gam)
    % Loop over number of users
    for ku = 1:length(U)
        % Loop over block sizes
        for kl = 1:length(L)

            % Evaluate Eq. (39)
            DeltaPout(ku,kl,kg) = DeltaOutage(gam(kg), mu2, U(ku),...
                L(kl), 'Quadrature', 30);
        end
    end
end

%-------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------
figure(1)

for ku = 1:length(U)
    plot(L, DeltaPout(ku,:,1), 'k-o','LineWidth',2, 'MarkerIndices',1:10:length(L))
    grid on; hold on
    plot(L, DeltaPout(ku,:,2), 'r--^','LineWidth',2, 'MarkerIndices',1:10:length(L))
end
l = legend("$\gamma = 0.5$", "$\gamma = 1$");
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('$L$', 'FontSize', 20, 'Interpreter','latex');
ylabel('$\Delta P_{out}$', 'FontSize', 20, 'Interpreter','latex');
ylim([0 0.3])
