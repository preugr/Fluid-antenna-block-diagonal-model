%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates and plots the OP achieved by slow-FAMA for 2D fluid antennas
% in terms of the SIR threshold. Clarke's correlation model is compared
% with the proposed block-diagonal approximation in Lemma 2 and the further
% approximation in Corollary 2. 
% MAKES FIGURE 13 IN THE PAPER
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
Nd = 20;                 % Port density per wavelength.

% Antenna size (wavelength normalized). Each row corresponds to a different
% fluid antenna (size along x-axis x size along z-axis)
W =  [2 1;
      4 2;
      5 3];             

U = [3 4 5 6 7 8];        % Number of users

Nsamples = 1e5;             % Number of samples for Monte-Carlo simulation

gam = 2;                    % SIR threshold (linear scale)

%-------------------------------------------------------------------------
% Computing OP
%-------------------------------------------------------------------------
% Clarke's correlation
pout_clarke = zeros(size(W,1), length(U));

% Block diagonal approximation OP in Lemma 2
pout_blocks = zeros(size(W,1),length(U));

% Approximated block diagonal OP in Corollary 2
pout_approx = zeros(size(W,1),length(U));

% Loop over fluid antenna area
for kw = 1:size(W,1)

    % User's feedback
    disp(['Outer loop iter ' num2str(kw) ' out of ' num2str(size(W,1))]);

    %-----------------------------------------------------------------
    % 3D Clarke's correlation
    %----------------------------------------------------------------
    N = W(kw,:) * Nd;
 
    % Antenna ports coordinates
    d = W(kw,:)./(N-1);   
    x_pos = (0:N(1,2)-1)*d(2);
    y_pos = (0:N(1,1)-1)*d(1);

    [X, Y] = ndgrid(x_pos, y_pos);
    xy_pos = [X(:) Y(:)];

    % Clarke's correlation matrix
    Sigma_clarke = SigmaIso3D(xy_pos);

    rho = eig(Sigma_clarke);
    rho = sort(rho,'descend');

    %----------------------------------------------------------------
    % Block diagonal correlation matrix approximation
    %-----------------------------------------------------------------
    Ntotal = N(1,1)*N(1,2);
    mu2 = 0.95;
    Num_eig = sum(rho > 1);
    % Algorithm 1
    L = BlockCorrelation(Ntotal, rho, Num_eig, mu2);

    % Loop over number of users
    for ku = 1:length(U)

        % User's feedback
        disp(['Inner loop iter ' num2str(ku) ' out of ' num2str(length(U))]);

        %-----------------------------------------------------------------
        % Outage Probabilities calculation/simulation
        %----------------------------------------------------------------- 
        % Clarke's outage probability
        pout_clarke(kw,ku) = SimOutage(Nsamples, gam, Sigma_clarke, U(ku));
    
        % Block outage probability -> Quadrature expression in Eq. (31)
        pout_blocks(kw,ku) = CalcOutage(gam, L, mu2, U(ku), 'Quadrature', 30);
    
        % Approximated blocks outage probability -> Quadrature expression
        % in Eq. (37)
        pout_approx(kw,ku) = ApproxOutage(gam, mu2, U(ku), L, 'Quadrature', 30);
    end

end


%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for kw = size(W,1):-1:1
    semilogy(U, pout_clarke(kw,:), 'k-*', 'linewidth', 2);
    hold on; grid on;
    semilogy(U, pout_blocks(kw,:), 'b-o', 'linewidth', 2);
    semilogy(U, pout_approx(kw,:), 'b--^', 'linewidth', 2);
end

l = legend("3D Clarke's", "Blocks (31)", "Blocks (37)");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('$U$', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{out}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([1e-4, 1])
xlim([U(1) U(end)])



