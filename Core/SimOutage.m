%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = SimOutage(Nsamples,gamma, Sigma, U) 
%
% Simulates the outage probability achieved by FAMA under Rayleigh fading
% with correlation matrix Sigma. If U == 1, then a single
% user case is simulated where gamma represents the SNR threshold instead 
% of the SIR threshold. Parameters:
%
% - Nsamples: number of Monte-Carlo simulations
% - gamma: SIR threshold (can be a vector) or SNR threshold if U == 1
% - U: number of users (scalar)
% - Sigma: correlation matrix of size NxN, where N is the number of ports.
%
% - pout: output vector (same size as gamma) with the outage probabilities.
%         pout(k) = P(SIR < gamma(k))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pout] = SimOutage(Nsamples,gamma, Sigma, U)

    % Pre-allocating
    pout = zeros(size(gamma));

    % Single user flag
    singleUser = (U == 1);

    % For the sake of efficiency, only a few dominant eigenvectors and
    % eigenvalues can be computed, specially useful as N becomes large. If
    % desired, the full eigendecomposition can be computed always as
    % eig(Sigma).
    if numel(Sigma) > 1e5
        [V,Lambda] = eigs(Sigma,100);
    else
        [V,Lambda] = eig(Sigma);
    end

    % Sort eigenvalues in descending order
    [lambda, index] = sort(diag(Lambda),'descend');

    % Only the eigenvalues larger than a small tolerance are stored
    lambda = lambda(lambda>1e-5);

    % Re-arranges the corresponding eigenvectors
    V = V(:,index(1:length(lambda)));

    % Batch size to avoid RAM issues. Nsamples are generated in batches of
    % size batchsize
    batchsize = 1e3;

    % Generate variables by batches
    for batch = 1:batchsize:Nsamples

        % Single user
        if singleUser
            X = GenVariables_SingleUser(batchsize, V, lambda);
            sir = max(X,[],1);
        % Multi-user
        else
            % Generate channels from desired user (X) and interference (Y)
            [Y, X] = GenVariables(batchsize, U, V, lambda);
    
            % Computes the SIR according to Eq. (24)
            sir = max(X./Y,[], 1);
        end
        % Computes outage probability
        for kg = 1:length(pout)
            pout(kg) = pout(kg) + sum(sir < gamma(kg));
        end
    end
    pout = pout/Nsamples;


    %---------------------------------------------------------------------
    %   Auxiliary function that generates correlated samples for the
    %   channel from/to the desired user and the interference. The channel
    %   is Rayleigh distributed with correlation matrix given by
    %   Sigma = V'*diag(lambda)*V
    %---------------------------------------------------------------------
    function [Y, X] = GenVariables(Nsamples, U, V, lambda)
    
        N = size(V,1);
    
        % Pre-allocating
        X = zeros(N, Nsamples);
        Y = zeros(N, Nsamples);
    
        % Square-root of correlation matrix
        L = V*diag(sqrt(lambda));

        % Number of independent Gaussians (rank of Sigma)
        Ngen = length(lambda);
    
        % Loop to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Correlated samples of channel from/to desired user
            X(:,junk_index) = abs(L*(randn(Ngen, length(junk_index)) + ...
                1i*randn(Ngen, length(junk_index)))/sqrt(2)).^2;
    
            % Correlated samples of channel from/to interferers
            for u = 1:(U-1)
                aux = abs(L*(randn(Ngen, length(junk_index)) + ...
                    1i*randn(Ngen, length(junk_index)))/sqrt(2)).^2;
                    Y(:,junk_index) = Y(:,junk_index) + aux; 
            end
    
        end
    
    end

    %---------------------------------------------------------------------
    %   Auxiliary function that generates correlated samples for the
    %   channel from/to the desired user. The channel
    %   is Rayleigh distributed with correlation matrix given by
    %   Sigma = V'*diag(lambda)*V
    %---------------------------------------------------------------------
    function [X] = GenVariables_SingleUser(Nsamples, V, lambda)
    
        N = size(V,1);
    
        X = zeros(N, Nsamples);
    
        L = V*diag(sqrt(lambda));

        Ngen = length(lambda);
    
        % Loop over junks to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            X(:,junk_index) = abs(L*(randn(Ngen, length(junk_index)) + ...
                1i*randn(Ngen, length(junk_index)))).^2;
    
        end
    
    end

end

