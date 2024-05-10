%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates and plots the approximated and real eigenvalues for both 
% Jakes's and Clarke's correlation models.
% MAKES FIG 5 IN THE PAPER
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
N = 120;                                    % Number of ports
W = [3 5];                                  % Lengths of fluid antenna (1D)

mu2 = 0.97;                                 % mu^2 for block approximation
rho_th = 1;                                 % Threshold for eigenvalues

colors = ['k', 'b', 'r', 'g', 'y', 'c'];    % Colours palette

%-------------------------------------------------------------------------
% Loop over W to compute eigenvalues (both real and approximated)
%-------------------------------------------------------------------------
for kw = 1:length(W)

    % Jakes's correlation matrix
    Sigma_jakes = toeplitz(besselj(0,2*pi*(0:N-1)*W(kw)/(N-1)));

    % Jakes's eigenvalues (sorted)
    rho_jakes = eig(Sigma_jakes);
    rho_jakes = sort(rho_jakes,'descend');

    % Clarke's (3D Isotropic) correlation matrix
    Sigma_clarke =  SigmaIso3D([0:W(kw)/(N-1):W(kw); zeros(1,N)]');

    % Clarke's eigenvalues (sorted)
    rho_clarke = eig(Sigma_clarke);
    rho_clarke = sort(rho_clarke,'descend');
    
    % Algorithm 1 (Block-diagonal approximation) for Jakes's model
    Num_eig = sum(rho_jakes > rho_th);
    L_jakes = BlockCorrelation(N, rho_jakes, Num_eig, mu2);

    % Algorithm 1 (Block-diagonal approximation) for Clarke's model
    Num_eig = sum(rho_clarke > rho_th);
    L_clarke = BlockCorrelation(N, rho_clarke, Num_eig, mu2);

    % Block covariance matrices --- NOTE: the eigenvalues can be calculated
    % directly from Eq. (21) using L and mu2. However, we explicitly build
    % the block-diagonal correlation matrix to cross-validate the result
    
    % Jakes approximation
    Sigma_blocks = zeros(N); 
    last_index = 0;
    for ksigma = 1:length(L_jakes)
        index = last_index + (1:L_jakes(ksigma));
        Sigma_blocks(index, index) = eye(L_jakes(ksigma)) + ...
                mu2*(1-eye(L_jakes(ksigma)));
        last_index = index(end);
    end

    rho_jakes_hat = eig(Sigma_blocks);
    rho_jakes_hat = sort(rho_jakes_hat,'descend');

    % Clarke approximation
    Sigma_blocks = zeros(N); 
    last_index = 0;
    for ksigma = 1:length(L_clarke)
        index = last_index + (1:L_clarke(ksigma));
        Sigma_blocks(index, index) = eye(L_clarke(ksigma)) + ...
                mu2*(1-eye(L_clarke(ksigma)));
        last_index = index(end);
    end

    rho_clarke_hat = eig(Sigma_blocks);
    rho_clarke_hat = sort(rho_clarke_hat,'descend');

    % Plotting
    figure(1) % Jakes's model
    t = 1:15;  % number of eigenvalues to plot
    plot(t, rho_jakes(t), colors(kw), 'linewidth',2,'DisplayName',...
        ['Jakes model W=' num2str(W(kw))]);
    hold on; grid on;
    plot(t, rho_jakes_hat(t), [colors(kw+2) '--'], 'linewidth',2,'DisplayName',...
        ['Block model W=' num2str(W(kw))])

    figure(2) % Clarke's model
    plot(t, rho_clarke(t), colors(kw), 'linewidth',2,'DisplayName',...
        ['Clarke model W=' num2str(W(kw))]);
    hold on; grid on;
    plot(t, rho_clarke_hat(t), [colors(kw+2) '--'], 'linewidth',2,'DisplayName',...
        ['Block model W=' num2str(W(kw))])


end

% Formatting the figures
figure(1)
l = legend();
set(gca, 'TickLabelInterpreter', 'latex','FontSize',18) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('n', 'FontSize', 18, 'Interpreter','latex');
ylabel('eigenvalues, $\rho_n$ and $\widehat{\rho}_n$', 'FontSize', 18, 'Interpreter','latex');

figure(2)
l = legend();
set(gca, 'TickLabelInterpreter', 'latex','FontSize',18) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('n', 'FontSize', 18, 'Interpreter','latex');
ylabel('eigenvalues, $\rho_n$ and $\widehat{\rho}_n$', 'FontSize', 18, 'Interpreter','latex');








