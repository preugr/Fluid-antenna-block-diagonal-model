%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = SimOutage_Blocks(Nsamples,gamma, U, mu, L)
%
% Simulates the outage probability achieved by FAMA under Rayleigh fading
% with block-diagonal correlation matrix as in Eq. (18). Parameters:
%
% - Nsamples: number of Monte-Carlo simulations
% - gamma: SIR threshold (can be a vector)
% - U: number of users (scalar)
% - mu: correlation factor (scalar) within each block (See Section III-B).
% - L: vector containing the size of each block (See Section III-B)
%
% - pout: output vector (same size as gamma) with the outage probabilities.
%         pout(k) = P(SIR < gamma(k)) according to Eq. (26)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pout] = SimOutage_Blocks(Nsamples,gamma, U, mu, L)

    % Pre-allocating
    pout = zeros(size(gamma));

    % Single user flag
    singleUser = (U == 1);

    % Batch size to avoid RAM issues. Nsamples are generated in batches of
    % size batchsize
    batchsize = 1e3;

    % Generate variables by batches
    for batch = 1:batchsize:Nsamples

        % Single user
        if singleUser
            X = GenVariables_SingleUser(batchsize, mu, L);
            sir = max((1-mu^2)*X,[],1);
        % Multiuser case
        else
            % Generate variables X and Y according to Eq. (26)
            [Y, X] = GenVariables(batchsize, U, mu, L);

            % Computes the SIR according to Eq. (26)
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
    function [Y, X] = GenVariables(Nsamples, U, mu, L)
    
        % Number of ports
        N = sum(L);
    
        % Pre-allocating
        X = zeros(N, Nsamples);
        Y = zeros(N, Nsamples);
    
        % Loop to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Loop over independent blocks
            for b = 1:length(L)
                b_index = sum(L(1:(b-1)))+1:sum(L(1:b));
                % Correlated samples of channel from/to desired user
                X(b_index,junk_index) = (randn(length(b_index),junk_size)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size)).^2 + ...
                        (randn(length(b_index),junk_size)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size)).^2;
    
                % Correlated samples of channel from/to interferers
                Y(b_index,junk_index) = sum((randn(length(b_index),junk_size,U-1)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size,U-1)).^2 + ...
                        (randn(length(b_index),junk_size,U-1)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size,U-1)).^2, 3);
            end
    
        end
    
    end

    function [X] = GenVariables_SingleUser(Nsamples, mu, L)
    
        % Number of ports
        N = sum(L);
    
        % Pre-allocating
        X = zeros(N, Nsamples);
    
        % Loop to avoid RAM issues
        junk_size = 1e3;
        for ks = 1:ceil(Nsamples/junk_size)
            junk_index = ((ks-1)*junk_size+1):min(junk_size*ks, Nsamples);
            % Loop over independent blocks
            for b = 1:length(L)
                b_index = sum(L(1:(b-1)))+1:sum(L(1:b));
                % Correlated samples of channel from/to desired user
                X(b_index,junk_index) = (randn(length(b_index),junk_size)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size)).^2 + ...
                        (randn(length(b_index),junk_size)...
                    + mu/sqrt(1-mu^2)*randn(1, junk_size)).^2;
    
            end
    
        end
    
    end

end

