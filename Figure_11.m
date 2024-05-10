%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates and plots the OP of slow-FAMA for 1D and 2D fluid antennas as
% their sizes vary. Clarke's correlation model is compared with the
% proposed approximation in Lemma 2 and the further approximation in
% Corollary 2.
% MAKES FIG 11 IN THE PAPER
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
% Number of ports (each row is a different fluid antenna)
N = [60 1;
     60 15;
     60 45];   

% Antenna size (wavelength normalized). Each row is a different fluid
% antenna
W = [4 0;
     4 1;
     4 3];                  

U = 6;                   % Number of users

Nsamples = 1e6;          % Number of samples for Monte-Carlo simulation

gamdB = linspace(-10,10,20);     % SIR threshold (dB)
gam = 10.^(gamdB/10);            % SIR threshold (linear scale)

%-------------------------------------------------------------------------
% Computing OP
%-------------------------------------------------------------------------
% Pre-allocating 

% Clarke's correlation
pout_clarke = zeros(size(W,1),length(gam));

% Block-diagonal approximation in Lemma 2
pout_blocks = zeros(size(W,1),length(gam));

% Further approximation in corollary 2
pout_blocks_approx = zeros(size(W,1),length(gam));

% Loop over fluid antenna area 
for kw = 1:size(W,1)

    % User feedback
    disp(['Iter ' num2str(kw) ' out of ' num2str(size(W,1))]);

    %-----------------------------------------------------------------
    % Clarke model
    %-----------------------------------------------------------------
    % Antenna ports coordinates
    if (N(kw,2) == 1)   % 1D case
        y_pos = (0:N(kw,1)-1)*W(kw,1)/(N(kw,1)-1);
        x_pos = 0;
    else                % 2D case 
        d = W(kw,:)./(N(kw,:)-1);   
        x_pos = (0:N(kw,2)-1)*d(2);
        y_pos = (0:N(kw,1)-1)*d(1);
    end
    
    [X, Y] = ndgrid(x_pos, y_pos);
    xy_pos = [X(:) Y(:)];

    % Clarke's correlation matrix
    Sigma_clarke = SigmaIso3D(xy_pos);

    rho = sort(eig(Sigma_clarke),'descend');
    
    %----------------------------------------------------------------
    % Block diagonal correlation matrix approximation
    %-----------------------------------------------------------------
    Ntotal = N(kw,1)*N(kw,2);
    mu2 = 0.97;
    Num_eig = sum(rho > 1);

    % Algorithm 1
    L = BlockCorrelation(Ntotal, rho, Num_eig, mu2);

    %-----------------------------------------------------------------
    % Outage Probabilities calculation/simulation
    %----------------------------------------------------------------- 
    % Clarke's model
    pout_clarke(kw,:) = SimOutage(Nsamples, gam, Sigma_clarke, U);

    % Block outage probability -> Evaluated through quadrature expression
    % in Eq. (31)
    pout_blocks(kw,:) = CalcOutage(gam, L, mu2, U, 'Quadrature', 30);

    % Approximated blocks outage probability -> Evaluated through
    % quadrature expression in Eq. (37)
    pout_blocks_approx(kw,:) = ApproxOutage(gam, mu2, U, L, 'Quadrature', 30);
end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for kw = 1:size(W,1)

    semilogy(gamdB, pout_clarke(kw,:), 'k-*', 'linewidth', 2);
    hold on; grid on;
    semilogy(gamdB, pout_blocks(kw,:), 'b-o', 'linewidth', 2);
    semilogy(gamdB, pout_blocks_approx(kw,:), 'b--^', 'linewidth', 2);
end

l = legend("3D Clarke's", "Blocks Eq. (29)", "Blocks Eq. (35)");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('$\gamma$ (dB)', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{out}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-5, 1])

