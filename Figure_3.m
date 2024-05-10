%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates and plots the number of eigenvalues surpassing a predefined
% threshold for both Jakes's and Clarke's correlation models in a 1D fluid 
% antenna. 
% MAKES FIGURE 3 IN THE PAPER
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
W = [1 2 3 4 5 6 7 8];                  % Normalized fluid antenna length
Ndensity = 20;                          % #ports per wavelength

%-------------------------------------------------------------------------
% Compute relevant eigenvalues, i.e., eigenvalues surpassing a threshold
%-------------------------------------------------------------------------
% Pre-allocating
relevant_jakes = zeros(length(W),1);
relevant_clarke = zeros(length(W),1);

% Loop over the different fluid antenna lengths
for kw = 1:length(W)

    % Total number of ports
    N = Ndensity * W(kw);

    % Ports positions (the fluid antenna is assumed to be along the y-axis)
    y_pos = (0:N-1)*W(kw)/(N-1);
    x_pos = 0;
    [X, Y] = ndgrid(x_pos, y_pos);
    xy_pos = [X(:) Y(:)];

    % Jakes's correlation matrix (Eq. (8))
    Sigma_jakes = toeplitz(besselj(0, 2*pi*(0:N-1)*W(kw)/(N-1)));

    % Clarke's correlation matrix (Eq. (7))
    Sigma_clarke = SigmaIso3D(xy_pos);

    % Threshold set so that dominant eigenvalues contain 99% of the power
    eig_th = N/100;    

    % Compute eigenvalues
    rho_jakes = eig(Sigma_jakes);
    rho_clarke = eig(Sigma_clarke);

    % Count dominant eigenvalues
    relevant_jakes(kw) = sum(rho_jakes > eig_th);
    relevant_clarke(kw) = sum(rho_clarke > eig_th);

end

%-------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------
figure(1)
plot(W, relevant_jakes, 'k-o', 'linewidth', 2, 'DisplayName', "Jakes's Model");
hold on; grid on;
plot(W, relevant_clarke, 'rv--', 'linewidth', 2, 'DisplayName', "Clarke's Model");

l = legend;

set(gca, 'TickLabelInterpreter', 'latex','FontSize',18) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('W (wavelength normalized)', 'FontSize', 18, 'Interpreter','latex');
ylabel('\# relevant eigenvalues', 'FontSize', 18, 'Interpreter','latex');

