%% Computation of illumination by a set of light sources
% 
% The set of light sources is defined in a structure array with fields:
% - HTM: the homogenous transformation matrix, defining position and
%        orientation;
% - 



%% Define receivers and compute values

% - Ar (receiver area) 
Ar = dx*dy;

rec_power = zeros(numel(yvals),numel(xvals));

for ix = 1:numel(xvals)
    for iy = 1:numel(xvals)
        % Compute receiver position
        HTM_R = Trans3(xvals(ix),yvals(iy),0);
        
        rp = 0;
        % for all emitters
        for e = Emitters
            % Compute received power
            rp = rp + e.Pt*H0_ER(e.HTM, HTM_R, e.m, Ar);
        end
        
        % Store received power
        % Note: the triples for surf(X,Y,Z) are defined as 
        % ( X(j), Y(i), Z(i,j) )
        % see doc surf
        rec_power(iy,ix) = rp;
    end
end


