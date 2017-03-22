

addpath('../ProjGeom/');

%% Define base values

% Room dimensions
Lx = 3;
Ly = 3;
Lz = 2; 


% Floor grid
% dx = 0.005;
% dy = 0.005; 

xvals = 0:dx:Lx;
yvals = 0:dy:Ly;

nx = numel(xvals);
ny = numel(yvals);



% Emitter position and orientation: 
HTM_E1 = Trans3(Lx/6, Ly/6, Lz);
HTM_E1 = HTM_E1 * RotX3(pi);
% Tilt the emitters
HTM_E1 = HTM_E1 * RotX3(pi/12) * RotY3(pi/12);

Emitters(1).HTM = HTM_E1;
Emitters(2).HTM = Trans3(0,2*Ly/3,0) * HTM_E1;
Emitters(3).HTM = Trans3(2*Lx/3,0,0) * HTM_E1;
Emitters(4).HTM = Trans3(2*Lx/3,2*Ly/3,0) * HTM_E1;

% Define default values
for i = 1:numel(Emitters)
    % Transmitted power
    Emitters(i).Pt = 1;
    % - m (Lambertian mode number)
    Emitters(i).m = 100;
end