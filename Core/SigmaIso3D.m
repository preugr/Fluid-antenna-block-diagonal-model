%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Sigma] = SigmaIso3D(xy_pos)
%
% Computes 3D Clarke's correlation matrix according to Eq. (11).
% Parameters:
%
% - xy_pos: rectangular coordinates of the ports composing the fluid 
%           antenna. The coordinates are normalized by the wavelength, and
%           can represent either a linear or a planar surface. Format is
%           xy_pos = [x_1 y_1;
%                     x_2 y_2;
%                       ...
%                     x_N y_N]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sigma] = SigmaIso3D(xy_pos)

    % Total number of ports
    Ntotal = size(xy_pos,1);

    Sigma = ones(Ntotal);
    for krow = 1:Ntotal
        for kcol = 1:Ntotal
            dist = sqrt(sum((xy_pos(kcol,:) - xy_pos(krow, :)).^2));            
            Sigma(krow, kcol) = sinc(2*dist);
        end
    end
end