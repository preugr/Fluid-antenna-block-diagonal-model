%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates and plots the number of eigenvalues surpassing a predefined
% threshold for Clarke's correlation model as well as for different
% correlation matrices where the azimuth angle domain is restricted. A 2D
% fluid antenna is assumed. 
% MAKES FIGURE 4 IN THE PAPER
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
% Without loss of generality, the fluid antenna is assumed to be located
% within the xz-plane
Wz = 4;                         % Dimension of fluid-antenna along z-axis
Wx = [0.5 1 2 3];               % Dimension of fluid-antenna along x-axis
Ndensity = 10;                  % Number of ports per wavelength
Nz = Ndensity*Wz;               % Total number of ports in the z direction

%-------------------------------------------------------------------------
% Compute relevant eigenvalues, i.e., eigenvalues surpassing a threshold
%-------------------------------------------------------------------------
% Pre-allocating
% Clarke's model
relevant_eig = zeros(length(Wx),1);               

% Azimuth span restricted to [pi/2 - phi_s/2, phi/2+phi_s/2]
relevant_eig_3pi4 = zeros(length(Wx),1);        % phi_s = 3pi/4
relevant_eig_pi2 = zeros(length(Wx),1);         % phi_s = pi/2
relevant_eig_pi4 = zeros(length(Wx),1);         % phi_s = pi/4

% Loop over fluid antenna area
for kw = 1:length(Wx)

    % User feedback
    disp(['Iter ' num2str(kw) ' out of ' num2str(length(Wx))]);

    % Total number of ports along x direction
    Nx = Ndensity * Wx(kw);

    % Ports coordinates
    z_pos = (0:Nz-1)*Wz/(Nz-1);
    x_pos = (0:Nx-1)*Wx(kw)/(Nx-1);
    [X, Z] = ndgrid(x_pos, z_pos);
    xz_pos = [X(:) Z(:)];

    % Clarke's correlation model
    Sigma_clarke = SigmaIso3D(xz_pos);

    % Different correlation matrices when azimuth span is restricted
    Sigma_3pi4= LimitedCorrelation_Iso(xz_pos, 3*pi/4, 30);
    Sigma_pi2= LimitedCorrelation_Iso(xz_pos, pi/2, 30);
    Sigma_pi4= LimitedCorrelation_Iso(xz_pos, pi/4, 30);

    % Threshold to determine dominant eigenvalues
    eig_th = 1;    

    % Computing eigenvalues
    rho_clarke = eig(Sigma_clarke);
    rho_3pi4 = eig(Sigma_3pi4);
    rho_pi2 = eig(Sigma_pi2);
    rho_pi4 = eig(Sigma_pi4);

    % Counting dominant eigenvalues
    relevant_eig(kw) = sum(rho_clarke > eig_th);
    relevant_eig_pi2(kw) = sum(rho_pi2 > eig_th);
    relevant_eig_3pi4(kw) = sum(rho_3pi4 > eig_th);
    relevant_eig_pi4(kw) = sum(rho_pi4 > eig_th);

end

%-------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------
figure(1)
plot(Wx*Wz, relevant_eig, 'k-o', 'linewidth', 2, 'DisplayName', "Clarke's ($\phi_{s} = 2\pi$)");
hold on; grid on;
plot(Wx*Wz, relevant_eig_3pi4, 'b-d', 'linewidth', 2, 'DisplayName', "$\phi_{s} = 3\pi/4$");
plot(Wx*Wz, relevant_eig_pi2, 'r-^', 'linewidth', 2, 'DisplayName', "$\phi_{s} = \pi/2$");
plot(Wx*Wz, relevant_eig_pi4, '-s', 'Color', [0.4660 0.6740 0.1880], 'linewidth', 2, 'DisplayName', "$\phi_{s} = \pi/4$");
l = legend;

set(gca, 'TickLabelInterpreter', 'latex','FontSize',18) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('$W_z \times W_x$', 'FontSize', 18, 'Interpreter','latex');
ylabel('\# relevant eigenvalues', 'FontSize', 18, 'Interpreter','latex');

%-------------------------------------------------------------------------
% Auxiliary function. Computes correlation matrix between antenna ports at
% coordinates indicated by xz_pos according to Eq. (5). The radiation
% pattern of the antennas is assumed to be isotropic (G(phi, theta) = 1),
% and the angle distribution is given by 
%     f_theta(theta) = sin(theta)/2         (Eq. (6))
%     phi (azimuth) is uniformly distributed between [pi/2 - phi_s/2, phi/2+phi_s/2]
%
% For the sake of efficiency, the double integral in Eq. (5) is computed
% through a nested application of the trapezoidal rule with K^2 total
% points.
% Parameters:
% - xz_pos: coordinates of the antenna ports in format
%           xz_pos = [x_1 z_1;
%                     x_2 z_2;
%                     ...
%                     x_N z_N];
% - phi_max: value of phi_s
% - K: points for each trapezoidal rule application
%-------------------------------------------------------------------------
function [Sigma] = LimitedCorrelation_Iso(xz_pos, phi_max, K)

      Ntotal = size(xz_pos,1);
      Sigma = zeros(Ntotal);
      k = 1:K-1;
      for nrow = 1:Ntotal
          for ncol = 1:Ntotal
              dist = (xz_pos(ncol,:) - xz_pos(nrow, :));  
              fun = @(theta, phi) exp(1i*2*pi*(sin(theta).*cos(phi)*dist(1)...
                  + cos(theta)*dist(2))).*sin(theta);
              
              Sigma(nrow,ncol) = (0.5*sum(fun(k*pi/K, pi/2-phi_max/2)) ...
                  + sum(fun(k*pi/K, pi/2-phi_max/2 + k'*phi_max/K),'all') ...
                  + 0.5*sum(fun(k*pi/K, pi/2+phi_max/2))) * phi_max*pi/(K^2*2*phi_max);
          end
      end

end





