%% Create receivers struct array

%% Initial Setup

addpath('../ProjGeom/');

% Define base values

% Room dimensions
Lx = 3;
Ly = 3;
Lz = 2; 


% Floor grid
dx = 0.05;
dy = 0.05; 

xvals = 0:dx:Lx;
yvals = 0:dy:Ly;

nx = numel(xvals);
ny = numel(yvals);

%% Create receivers

% Create the receiver structure:
Receiver_t = struct('HTM',{},'Ar',{},'Ts',{},'n',{},'Psi',{},'Pr',{});

% Create one element with default values
Receiver_t(1).Ar = dx*dy;
Receiver_t(1).Ts = 1;
Receiver_t(1).n = 1;
Receiver_t(1).Psi = pi/2;
Receiver_t(1).Pr = 0;

% Replicate to create the receivers array:
Receivers = repmat(Receiver_t,1,10);

HTM_R_base = Trans3(Lx/6, Ly/6, 0);

Receivers(1).HTM = HTM_R_base;
Receivers(2).HTM = Trans3(0,2*Ly/3,0) * HTM_R_base;
Receivers(3).HTM = Trans3(2*Lx/3,0,0) * HTM_R_base;
Receivers(4).HTM = Trans3(2*Lx/3,2*Ly/3,0) * HTM_R_base;

Receivers(5).HTM = Trans3(0, 0, 2);
Receivers(6).HTM = Trans3(0, Ly, 2);
Receivers(7).HTM = Trans3(Lx, 0, 2);
Receivers(8).HTM = Trans3(Lx, Ly, 2);

Receivers(9).HTM = Trans3(Lx/2, Ly/2, 0);
Receivers(10).HTM = Trans3(Lx/2, Ly/2, 2);

