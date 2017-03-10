
addpath('../ProjGeom/');

% Define base values

% Room dimensions
Lx = 3;
Ly = 3;
Lz = 2; 

dx = .05; 
dy = .05;

% Emitter position and orientation: 
HTM_E1 = Trans3(Lx/2, Ly/4, Lz);
HTM_E1 = HTM_E1 * RotX3(pi+pi/12) * RotY3(pi/12);

HTM_E2 = Trans3(0,Ly/2,0) * HTM_E1;

% - m (Lambertian mode number)
% - Ar (receiver area) 
m = 100;
Ar = 1e-4;


xvals = 0:dx:Lx;
yvals = 0:dy:Ly;

rec_power = zeros(numel(xvals),numel(yvals));

% Compute values 
for ix = 1:numel(xvals)
    for iy = 1:numel(xvals)
        % Compute receiver position
        HTM_R = Trans3(xvals(ix),yvals(iy),0);
        
        % Compute received power
        rp1 = H0_ER(HTM_E1, HTM_R, m, Ar);
        rp2 = H0_ER(HTM_E2, HTM_R, m, Ar);
        
        % Store received power
        % Note: the triples for surf(X,Y,Z) are defined as 
        % ( X(j), Y(i), Z(i,j) )
        % see doc surf
        rec_power(iy,ix) = rp1 + rp2;
    end
end

% Show the results
figure
plot3Drefaxis(HTM_E1)
axis equal
hold
plot3Drefaxis(HTM_E2)

h = surf(xvals,yvals,rec_power);
h.MeshStyle = 'none';
% view(2);
% figure(gcf);
axis equal

