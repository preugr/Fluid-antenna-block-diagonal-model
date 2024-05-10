%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [L] = BlockCorrelation(N, rho, Num_eig, mu2)
%
% Implements Algorithm 1 (See Section III_B), which calculates the 
% different block sizes L_b to approximate the set of dominant eigenvalues 
% indicated in rho. 
%
% Parameters:
% - N: number of ports in fluid antenna (scalar)
% - rho: vector containing the set of eigenvalues of the target correlation
%        matrix
% - Num_eig: number of eigenvalues (dominant) in rho that want to be
%            approximated (scalar). 
% - mu2: mu^2 (scalar)
%
% - L: vector (size 1xNum_eig) containing the block sizes L_b 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L] = BlockCorrelation(N, rho, Num_eig, mu2)
    
    % Dominant set of eigenvalues
    lambda_eff = rho(1:Num_eig);
    
    %------------------------------ Algorithm 1 --------------------------
    L = zeros(Num_eig,1);
    counter = 0;
    index = ones(length(L),1);
    
    while counter < N && sum(index) > 0
        L(index==1) = L(index==1) + ones(sum(index),1);
        counter = sum(L);
    
        index(index==1) = abs((L(index==1)-1)*mu2 + 1 - lambda_eff(index==1)) >= ...
                abs((L(index==1))*mu2 + 1 - lambda_eff(index==1));
    
    end
end

