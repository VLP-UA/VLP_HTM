% Demonstrates the usage of HTM based VLP code

%% Prepare the workspace

% add the path to the projective geometry functions
addpath('../ProjGeom');

% Define the model default values

% Emitters
n_Emitters = 3;         % Number of emitters
Pt = 1 ;                % Transmitted power
m = 5 ;                 % Lambertian mode number

% Receivers:
Np = 10;                % Number of parallels in the sensor
Nm = 10;                % Number of meridians in the sensor
n_Receivers = Np*Nm;    % Number of receivers
Ar = 0.01;              % Active receiving area
Ts = 1;                 % Optical filter gain
n = 1;                  % Receiver's internal refractive index
Psi = pi/6;             % Hemi-Fov



%% Create and populate data structures

% Create the emitters array:
Emitters = newEmitters(n_Emitters,Pt, m);

% Create the receiver structure:
Receivers = newReceivers(n_Receivers,Ar, Ts, n, Psi);

%% Place emitters and receivers in 3D space

% Considering a room with 3mx3mx2m (WxLxH).

% Emitters are placed at the ceiling, looking down and in a circle of
% radius=0.5m

Em_Base_HTM = Trans3(1.5,1.5,2)*RotX3(pi);      % Base HTM at the center of the ceiling. 

Emitters(1).HTM = Em_Base_HTM*Trans3(0.5,0,0);
Emitters(2).HTM = Em_Base_HTM*RotZ3(2*pi/3)*Trans3(0.5,0,0);
Emitters(3).HTM = Em_Base_HTM*RotZ3(4*pi/3)*Trans3(0.5,0,0);

% Receivers are organized in Parallel and Meridians arragement of photo
% detectors, with Nm Meridians and 3 Parallels, in a sphere with
% radius=0.25
PDSensor = vlpCreateSensorParMer(Receivers, Np, Nm, 0.25);

% The sensor is place at point (X,Y) = (0.5,0)
PDSensor = vlpMoveSensor(PDSensor,Trans3(0.5,0,0));

%% Compute A, B and F matrices

A = Receiver2A(PDSensor);
B = Emitter2B(Emitters);
F = EmitRec2F(Emitters, PDSensor);

%% Compute the received power

Pt_v = [Emitters.Pt]';
Pr = A*F*B*Pt_v;

% Update the PDSensor received power
for i=1:numel(PDSensor)
    PDSensor(i).Pr = Pr(i);
end


%% Plot Emitters and receivers

% Plot the emitters and receivers position
figure;
PlotHTMArray(Emitters);
PlotHTMArray(PDSensor);
return
% Plot with received power intensity
figure
PlotHTMArray(Emitters);
PlotHTMArrayPr(PDSensor);
