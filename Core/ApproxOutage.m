%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pout] = ApproxOutage(gamma_v, rho, U, L, method, order)
%
% Calculates approximated OP of slow-FAMA for the block-diagonal correlation
% model as in Corollary 2. Parameters:
%
% - gamma_v: vector containing the SIR thresholds in linear scale.
% - rho: squared value of the correlation coefficient mu used in the 
%        block-diagonal approximation, i.e., rho = mu^2
% - U: number of users (scalar)
% - L: vector containing the block sizes of the correlation approximation
%      (See Section III-B)
% - method: string that specifies the theoretical expression used to
%           compute the OP. 
%           If method == 'Integral', then the single integral expression 
%           in Eq. (31) is used. 
%           If method == 'Quadrature', the integral is solved by quadrature
%           as in Eq. (37).
% - order: order of the quadrature approximation 
%          (only if method == 'Quadrature' )
%
% - pout: vector containing the approximated OP (same size as gamma_v), 
%         i.e., pout(k) = P(SIR < gamma_v(k)) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pout] = ApproxOutage(gamma_v, rho, U, L, method, order)

    % Number of blocks
    B = length(L);

    %---------------- Direct integration method (Eq. (35)) -------------
    if strcmp(method, 'Integral')
        % Pre-allocating
        pout = ones(1,length(gamma_v));

        for kmu = 1:B

            % Define integrand with particularized values
            fint = @(rtil_b) rtil_b.^(U-2).*exp(-(rtil_b + (sqrt(gamma_v.*rtil_b)...
                + ((U-1.5)*sqrt((gamma_v+1)*(1-rho))/sqrt(rho) - ...
                    (L(kmu)-1)*sqrt(gamma_v.*rtil_b/(2*pi)))./...
                    ((L(kmu)-1)*(U-1.5)/sqrt(2*pi) + ...
                    sqrt(rho*gamma_v.*rtil_b/(1-rho)./(1+gamma_v)))).^2)/2);
        
            pout = pout .* ...
                    (1-(integral(fint, 0, inf, 'ArrayValued',true))...
                        /(2^(U-1)*gamma(U-1)));

        end
    
    %---------------- Quadrature evaluation (Eq. (37)) -------------
    elseif strcmp(method, 'Quadrature')

        t = sym('t');
        % Computing roots of generalized Laguerre polynomials
        L_n_alpha = laguerreL(order, U-2, t);
        xtil = roots(sym2poly(L_n_alpha));
        wtil = real(exp(gammaln(order+U-1) +log(xtil) - gammaln(order+1) - ...
            2*log(order+1) - 2*log(laguerreL(order+1, U-2, xtil))));

        % Compute OP
        pout = zeros(length(L),length(gamma_v));
        for kmu = 1:B
            for ki = 1:order
                pout(kmu,:) = pout(kmu,:) + wtil(ki)* ...
                    exp(-((sqrt(gamma_v.*2*xtil(ki))...
                + ((U-1.5)*sqrt((gamma_v+1)*(1-rho))/sqrt(rho) - ...
                    (L(kmu)-1)*sqrt(gamma_v.*2*xtil(ki)/(2*pi)))./...
                    ((L(kmu)-1)*(U-1.5)/sqrt(2*pi) + ...
                    sqrt(rho*gamma_v.*2*xtil(ki)/(1-rho)./(1+gamma_v)))).^2)/2);
            end
        end

        pout = prod(1 - pout/gamma(U-1),1); 
       
    else
        error('Integration method must be indicated')
    end

end

