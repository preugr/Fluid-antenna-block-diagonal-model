%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = CalcOutage(gamma_v, L, rho, U, method, order)
%
% Calculates theoretical OP of slow-FAMA for the block-diagonal correlation
% approximation as in Eq. (29). If the number of users is set to 1, then a 
% single user fluid antenna system is considered, and the OP acccording 
% to Eq. (43) is computed, where gamma_v becomes the SNR (and not SIR)
% threshold. 
% 
% Parameters:
%
% - gamma_v: vector containing the SIR thresholds in linear scale. If U ==
%            1, then gamma_v represents the SNR threshold.
% - L: vector containing the block sizes of the correlation approximation
%      (See Section III-B)
% - rho: squared value of the correlation coefficient mu used in the 
%        block-diagonal approximation, i.e., rho = mu^2
% - U: number of users (scalar). If U ==1, then Eq. (43) is computed
% - method: string that specifies the theoretical expression used to
%           compute the OP. 
%           If method == 'Integral', then the double integral expression 
%           in Eq. (29) is used. 
%           If method == 'Quadrature', the integral is solved by quadrature
%           as in Eq. (31).
% - order: order of the quadrature approximation 
%          (only if method == 'Quadrature' )
%
% - pout: vector containing the OP (same size as gamma_v), i.e., 
%         pout(k) = P(SIR < gamma_v(k))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pout] = CalcOutage(gamma_v, L, rho, U, method, order)

    %--------------- Single user case -> Evaluation of Eq. (43) ----------
    %---------------------------------------------------------------------
    if (U == 1 || nargin < 4)
        % Pre-allocating
        pout = ones(1,length(gamma_v));

        for kb = 1:length(L)

        fun = @(rb) 0.5*exp(-rb/2).*...
            (1-marcumq(sqrt(rho*rb./(1-rho)),sqrt(gamma_v./(1-rho)),1)).^L(kb);

        pout = pout.*integral(fun, 0, inf, 'ArrayValued',true);

        end
        

        %--------------- Multi user case -------------------------------------
        %---------------------------------------------------------------------
    else
        %---------------- Direct integration method (Eq. (29)) -------------
        if strcmp(method, 'Integral')
            % Pre-allocating
            pout = ones(1,length(gamma_v));
    
            % Loop over SIR thresholds
            for kg = 1:length(pout)
        
                % Loop over blocks
                for kmu = 1:length(L)
        
                    % Define integrand with particularized values
                    fint = @(r_b, rtil_b) rtil_b.^(U-2).*exp(-(r_b+rtil_b)/2)...
                             .*Gfun(r_b, rtil_b, gamma_v(kg), U, rho).^L(kmu);
    
                    % Perform numerical double integration
                    pout(kg) = pout(kg) * ...
                            (integral2(fint, 0, inf, 0, inf))/(2^U*gamma(U-1));
        
                end
        
            end
        %---------------- Quadrature approximation (Eq. (31)) -------------
        elseif strcmp(method, 'Quadrature')
    
            % Pre-compute Laguerre polynomials roots
            t = sym('t');
            L_n = laguerreL(order, t);
            x = roots(sym2poly(L_n));
            w = x./((order+1)^2 * (laguerreL(order+1, x)).^2);
    
            % Pre-compute generalized Laguerre polynomials roots
            L_n_alpha = laguerreL(order, U-2, t);
            xtil = roots(sym2poly(L_n_alpha));
            wtil = real(exp(gammaln(order+U-1) +log(xtil) - gammaln(order+1) - ...
                2*log(order+1) - 2*log(laguerreL(order+1, U-2, xtil))));
            
            % Store each summand in a 3D tensor
            Gmatrix = zeros(order, order, length(gamma_v));
            for ki = 1:order
                for kj = 1:order
                    Gmatrix(ki,kj,:) = Gfun(2*x(kj), 2*xtil(ki), gamma_v, U, rho);
                end
            end
    
            % Computes outage probability
            pout = ones(length(gamma_v),1);
            for kmu = 1:length(L)
                pout = pout.* squeeze(sum(wtil.*(Gmatrix.^L(kmu)).*w.',[1 2]))/gamma(U-1);
            end
           
        else
            error('Integration method must be indicated')
        end
    end

    %---------------------------------------------------------------------
    % Auxiliary function to compute G as in Eq. (30)
    %---------------------------------------------------------------------
    function [G] = Gfun(r_b, rtil_b, gam, U, rho)

        summatory = 0;

        for k = 0:(U-2)
            for j = 0:(U-k-2)
                summatory = summatory + exp(gammaln(U-k-1) - gammaln(U-j-k-1)...
                    - gammaln(j+1) + ((j+k)/2)*log(r_b./rtil_b)) .* (gam+1).^k ...
                    .* gam.^((j-k)/2) .* ...
                    besseli(j+k, rho*sqrt(gam.*r_b.*rtil_b)/(1-rho)./(gam+1),1)...
                    .* exp((rho/(1-rho)./(gam+1)).*...
                        (sqrt(gam.*r_b.*rtil_b)-(gam.*rtil_b+r_b)/2));
            end
        end

        G = marcumq(sqrt(rho*gam.*rtil_b/(1-rho)./(gam+1)), ...
                       sqrt(rho*r_b/(1-rho)./(gam+1)),U-1) - ...
               (1./(gam+1)).^(U-1).*summatory;

        G(isnan(G)) = 0; 

    end


end

