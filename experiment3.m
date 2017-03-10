%% Computation of illumination by a set of light sources
% 
% The set of light sources is defined in a structure array with fields name
% (optional) and HTM (the homogenous transformation matrix)


addpath('../ProjGeom/');

% Define base values

% Room dimensions
Lx = 3;
Ly = 3;
Lz = 2; 

dx = .05; 
dy = .05;

% Emitter position and orientation: 
HTM_E1 = Trans3(Lx/6, Ly/6, Lz);
HTM_E1 = HTM_E1 * RotX3(pi);
% Tilt the emitters
%HTM_E1 = HTM_E1 * RotX3(pi/12) * RotY3(pi/12);

Emitters(1).HTM = HTM_E1;
Emitters(2).HTM = Trans3(0,2*Ly/3,0) * HTM_E1;
Emitters(3).HTM = Trans3(2*Lx/3,0,0) * HTM_E1;
Emitters(4).HTM = Trans3(2*Lx/3,2*Ly/3,0) * HTM_E1;

% - m (Lambertian mode number)
% - Ar (receiver area) 
m = 1;
Ar = 1e-4;


xvals = 0:dx:Lx;
yvals = 0:dy:Ly;

rec_power = zeros(numel(yvals),numel(xvals));

%% Compute values 
for ix = 1:numel(xvals)
    for iy = 1:numel(xvals)
        % Compute receiver position
        HTM_R = Trans3(xvals(ix),yvals(iy),0);
        
        rp = 0;
        % for all emitters
        for e = Emitters
            % Compute received power
            rp = rp + H0_ER(e.HTM, HTM_R, m, Ar);
        end
        
        % Store received power
        % Note: the triples for surf(X,Y,Z) are defined as 
        % ( X(j), Y(i), Z(i,j) )
        % see doc surf
        rec_power(iy,ix) = rp;
    end
end

%% Show the results

figure

% Plot the light sources HTM reference frames
for e = Emitters
    plot3Drefaxis(e.HTM)
    hold on
    axis equal
end

h = surf(xvals,yvals,rec_power);
h.MeshStyle = 'none';
% view(2);
% figure(gcf);
axis equal

