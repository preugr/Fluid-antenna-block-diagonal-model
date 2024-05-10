%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the evolution of OP as a 1D fluid antenna is saturated (i.e., more
% ports are added while keeping the aperture constant). Clarke's correlation
% model is compared with the proposed block-diagonal approximation.
% MAKES FIGURE 9 IN THE PAPER
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
Niid = [4 8 16 32];                % Number of antennas in i.i.d. case

W = [2 4 8];                       % Antenna size (wavelength normalized)
U = 3;                             % Number of users

Nsamples = 1e6;              % Number of samples for Monte-Carlo simulation

gam = 1;                     % SIR threshold (linear scale)

%-------------------------------------------------------------------------
% Computing OP
%-------------------------------------------------------------------------
% Pre-allocating
pout_clarke = zeros(length(W),length(Nd));
pout_blocks = zeros(length(W),length(Nd));

% Loop over port density
for ksize = 1:length(Nd)

    % Loop over fluid antenna length
    for kw = 1:length(W)

        % User feedback
        disp(['Iter ' num2str((ksize-1)*length(W)+kw) ' out of '...
            num2str(length(Nd)*length(W))]);

        %-----------------------------------------------------------------
        % 3D Clarke's correlation
        %----------------------------------------------------------------
        N = W(kw) * Nd(ksize);
     
        Sigma_clarke = toeplitz(sinc(2*(0:(N-1))*W(kw)/(N-1)));
    
        rho = sort(eig(Sigma_clarke),'descend');
    
        %----------------------------------------------------------------
        % Block diagonal correlation matrix approximation
        %-----------------------------------------------------------------
        mu2 = 0.97;
        Num_eig = sum(rho > N/100);
        L = BlockCorrelation(N, rho, Num_eig, mu2);
    
        %-----------------------------------------------------------------
        % Outage Probabilities calculation/simulation
        %----------------------------------------------------------------- 
        % Clarke outage probability
        pout_clarke(kw,ksize) = SimOutage(Nsamples, gam, Sigma_clarke, U);
    
        % Block outage probability
        pout_blocks(kw,ksize) = CalcOutage(gam, L, mu2, U, 'Quadrature', 30);
        
    end
end

%-----------------------------------------------------------------
% Independent antenna case -> Eq. (38)
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

for kw = length(W):-1:1
    semilogy(Nd, pout_clarke(kw,:), 'k-*', 'linewidth', 2);
    hold on; grid on;
    semilogy(Nd, pout_blocks(kw,:), 'b-o', 'linewidth', 2);
end

l = legend("", "", "", "", "3D Clarke's", "Blocks (29)");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 20, 'Interpreter','latex');
xlabel('$N/W$', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{out}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([8e-5, 1])
xlim([Nd(1) Nd(end)])




