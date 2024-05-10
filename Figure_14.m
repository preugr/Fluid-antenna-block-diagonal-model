%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates and plots the OP achieved by a single-user fluid antenna system
% in terms of the SNR threshold. Different fluid antenna sizes (either 1D 
% and 2D) are compared, and Clarke's correlation model is compared with the
% proposed block diagonal approximation
% MAKES FIGURE 14 IN THE PAPER
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
Nden = [2 4 6 8 10 12 14 16 18 20];     % # ports per wavelength

% Antenna size (wavelength normalized). Each row represents a different
% fluid antenna
W = [3 1;
     3 2;
     3 3];                  

Nsamples = 5e5;           % Number of samples for Monte-Carlo simulation
 
gam = 5;                  % SNR threshold (linear scale)

%-------------------------------------------------------------------------
% OP computation
%-------------------------------------------------------------------------
% Pre-allocating
pout_clarke = zeros(size(W,1),length(Nden));
pout_blocks = zeros(size(W,1),length(Nden));

% loop over fluid antenna areas
for kw = 1:size(W,1)

    % Loop over port density
    for kn = 1:length(Nden)

        % User feedback
        disp(['Iter ' num2str((kw-1)*length(Nden)+kn) ' out of '...
            num2str(size(W,1)*length(Nden))]);

        %-----------------------------------------------------------------
        % 3D Isotropic correlation
        %-----------------------------------------------------------------
        if (W(kw,2) == 1)               % 1D case
            N = W(kw,1)*Nden(kn);
            y_pos = (0:N-1)*W(kw,1)/(N-1);
            x_pos = 0;
        
            [X, Y] = ndgrid(x_pos, y_pos);
            xy_pos = [X(:) Y(:)];
        
            Sigma_clarke = SigmaIso3D(xy_pos);
        
            rho = sort(eig(Sigma_clarke),'descend');
        else                            % 2D case
            N = Nden(kn)*W(kw,:);
            d = W(kw,:)./(N-1);   
            x_pos = (0:N(2)-1)*d(2);
            y_pos = (0:N(1)-1)*d(1);
            
            [X, Y] = ndgrid(x_pos, y_pos);
            xy_pos = [X(:) Y(:)];
            
            Sigma_clarke = SigmaIso3D(xy_pos);
            
            rho = sort(eig(Sigma_clarke),'descend');

        end

        %----------------------------------------------------------------
        % Block diagonal correlation matrix approximation
        %-----------------------------------------------------------------
        Ntotal = prod(N);
        if W(kw,2) == 1
            mu2 = 0.95;
        else
            mu2 = 0.98;
        end
  
        Num_eig = sum(rho > 1);

        % Algorithm 1
        L = BlockCorrelation(Ntotal, rho, Num_eig, mu2);

        %-----------------------------------------------------------------
        % Outage Probabilities calculation/simulation
        %----------------------------------------------------------------- 
        % Clarke's model
        pout_clarke(kw,kn) = SimOutage(Nsamples, gam, Sigma_clarke, 1);

        % Block model
        pout_blocks(kw,kn) = CalcOutage(gam, L, mu2, 1);

    end
end

%---------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------- 
figure(1)

for kw = 1:size(W,1)

    semilogy(Nden, pout_clarke(kw,:), 'k-o', 'linewidth', 2);
    hold on; grid on;
    semilogy(Nden, pout_blocks(kw,:), 'b-o', 'linewidth', 2);
end

l = legend("Clarke's", "Blocks Eq. (43)");

set(gca, 'TickLabelInterpreter', 'latex','FontSize',20) 
set(l, 'FontSize', 18, 'Interpreter','latex');
xlabel('\# ports per wavelength', 'FontSize', 20, 'Interpreter','latex');
ylabel('$P_{out}(\gamma)$', 'FontSize', 20, 'Interpreter','latex');
ylim([5e-5, 1])
