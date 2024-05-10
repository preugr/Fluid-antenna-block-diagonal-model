%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DeltaPout] = DeltaOutage(gamma_v, rho, U, L, method, order)
%
% Evaluates the gain, in terms of OP, achieved by each block in the 
% block-diagonal correlation model according to Eq. (39). Parameters:
%
% - gamma_v: vector containing the SIR thresholds in linear scale.
% - rho: squared value of the correlation coefficient mu used in the 
%        block-diagonal approximation, i.e., rho = mu^2
% - U: number of users (scalar)
% - L: vector containing the block sizes of the correlation approximation
%      (See Section III-B)
% - method: string that specifies the theoretical expression used to
%           compute the OP. 
%           If method == 'Integral', then the double integral expression 
%           in Eq. (39) is used. 
%           If method == 'Quadrature', Eq. (39) is solved by quadrature.
% - order: order of the quadrature approximation 
%          (only if method == 'Quadrature' )
%
% - DeltaPout: vector containing the result of Eq. (39)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DeltaPout] = DeltaOutage(gamma_v, rho, U, L, method, order)

    %---------------------------------------------------------------------
    % INTEGRAL METHOD -> Direct evaluation of Eq. (39)
    %---------------------------------------------------------------------
    if strcmp(method, 'Integral')

        % Define integrand with particularized values
        fint = @(rtil_b) rtil_b.^(U-2).*exp(-(rtil_b + (sqrt(gamma_v.*rtil_b)...
            + ((U-1.5).*sqrt((gamma_v+1)*(1-rho))/sqrt(rho) - ...
                (L-1).*sqrt(gamma_v.*rtil_b/(2*pi)))./...
                ((L-1).*(U-1.5)/sqrt(2*pi) + ...
                sqrt(rho*gamma_v.*rtil_b/(1-rho)./(1+gamma_v)))).^2)/2);
    
        DeltaPout = (integral(fint, 0, inf, 'ArrayValued',true))...
                    /(2^(U-1)*gamma(U-1)) - 1./((gamma_v+1).^(U-1));
    
    %---------------------------------------------------------------------
    % QUADRATURE
    %---------------------------------------------------------------------
    elseif strcmp(method, 'Quadrature')

        % Calculate roots of modified Laguerre polynomials
        t = sym('t');

        L_n_alpha = laguerreL(order, U-2, t);
        xtil = roots(sym2poly(L_n_alpha));
        wtil = real(exp(gammaln(order+U-1) +log(xtil) - gammaln(order+1) - ...
            2*log(order+1) - 2*log(laguerreL(order+1, U-2, xtil))));

        % Compute Eq. (39) through quadrature
        DeltaPout = zeros(1,length(gamma_v));
        for ki = 1:order
            DeltaPout = DeltaPout + wtil(ki)* ...
                exp(-((sqrt(gamma_v.*2*xtil(ki))...
            + ((U-1.5)*sqrt((gamma_v+1)*(1-rho))/sqrt(rho) - ...
                (L-1)*sqrt(gamma_v.*2*xtil(ki)/(2*pi)))./...
                ((L-1)*(U-1.5)/sqrt(2*pi) + ...
                sqrt(rho*gamma_v.*2*xtil(ki)/(1-rho)./(1+gamma_v)))).^2)/2);
        end

        DeltaPout =  DeltaPout/gamma(U-1) - 1./((gamma_v+1).^(U-1)); 
       
    else
        error('Integration method must be indicated')
    end

end

