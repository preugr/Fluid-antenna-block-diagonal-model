%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates and plots the OP achieved by slow-FAMA as the 2D fluid antenna
% is saturated (i.e., the number of ports is increased while keeping 
% constant the aperture). Clarke's correlation model is compared with the
% proposed block-diagonal approximation.
% MAKES FIG 12 IN THE PAPER
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
Nd = [2 6 12 18 24 30];     % Port density per wavelength.

Niid = [64 128 256 512];    % Number of i.i.d. antennas as baseline

% Antenna size (wavelength normalized). Each row corresponds to a different
% fluid antenna (size along x-axis x size along z-axis)
W = [3 2;
     4 3;
     6 4];                       

U = 7;                             % Number of users

Nsamples = 5e5;                    % Number of samples for Monte-Carlo simulation

gam = 1;                           % SIR threshold (linear scale)

%-------------------------------------------------------------------------
% Computing OP
%-------------------------------------------------------------------------
% Pre-allocating

% Clarke's correlation model
pout_clarke = zeros(size(W,1),length(Nd));

% Block diagonal approximation OP in Lemma 2
pout_blocks = zeros(size(W,1),length(Nd));

% Loop over port density
for ksize = 1:length(Nd)

    % Loop over fluid antenna area
    for kw = 1:size(W,1)

        % User's feedback
        disp(['Iteration ' num2str((ksize-1)*size(W,1) + kw)...
            ' out of ' num2str(size(W,1)*length(Nd))]);

        %-----------------------------------------------------------------
        % 3D Clarke's correlation
        %----------------------------------------------------------------
        N = W(kw,:) * Nd(ksize);
     
        % Antenna ports coordinates
        d = W(kw,:)./(N-1);   
        x_pos = (0:N(1,2)-1)*d(2);
        y_pos = (0:N(1,1)-1)*d(1);
    
        [X, Y] = ndgrid(x_pos, y_pos);
        xy_pos = [X(:) Y(:)];

        % Clarke's correlation matrix
        Sigma_clarke = SigmaIso3D(xy_pos);

        rho = sort(eig(Sigma_clarke),'descend');
    
        %----------------------------------------------------------------
        % Block diagonal correlation matrix approximation
        %-----------------------------------------------------------------
        Ntotal = N(1,1)*N(1,2);
        mu2 = 0.96;
        Num_eig = sum(rho > 1);
        % Algorithm 1
        L = BlockCorrelation(Ntotal, rho, Num_eig, mu2);

        %-----------------------------------------------------------------
        % Outage Probabilities calculation/simulation
        %----------------------------------------------------------------- 
        % Clarke's outage probability
        pout_clarke(kw,ksize) = SimOutage(Nsamples, gam, Sigma_clarke, U);

        % Block outage probability -> Evaluated through quadrature
        % expression in Eq. (31)
        pout_blocks(kw,ksize) = CalcOutage(gam, L, mu2, U, 'Quadrature', 30);
        
    end
end

%-----------------------------------------------------------------
% Independent antenna case
%----------------------------------------------------------------- 
pout_iid = (1 - 1/(gam+1)^(U-1)).^Niid;

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for kiid = 1:length(Niid)
    semilogy(Nd, pout_iid(kiid)*ones(length(Nd),1), 'r-.', 'linewidth', 2);
    hold on;
end

for kw = size(W,1):-1:1
    semilogy(Nd, pout_clarke(kw,:), 'k-*', 'linewidth', 2);
    hold on; grid on;
    semilogy(Nd, pout_blocks(kw,:), 'b-o', 'linewidth', 2);
end

l = legend("", "", "", "", "3D Clarke's", "Blocks (31)");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('$\#$ ports $/$ wavelength', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{out}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-4, 1])
xlim([Nd(1) Nd(end)])


