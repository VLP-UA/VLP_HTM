%RUN_MULT_PARAMETERS_DAUGHTER- run a give set of parameters, refer to
%run_mult_parameters to run this script, params struct needs to be
%initialized first


%% Create the light emitters
broken = 0;
% Emitters
n_Emitters = params.n_Emitters;         % Number of emitters
Pb = params.Pb;                % Transmitted power
Ps = params.Ps;
m = params.m;                  % Lambertian mode number

Emitters = newEmitters(n_Emitters,Pb,Ps,m);

% Room dimensions
% Considering a room with 3mx3mx2m (WxLxH).
W = params.W;
L = params.L;
H = params.H;

offset = 5;

Em_Base_HTM = Trans3(0,0,H)*RotX3(pi);      % Base HTM at the center of the ceiling.

tempx = linspace(offset,W-offset,sqrt(n_Emitters));
tempy = linspace(offset,L-offset,sqrt(n_Emitters));
[temp_x,temp_y] = meshgrid(tempx,tempy);

temp_x=temp_x(:);
temp_y=temp_y(:);
% Light emitters placed at ceiling, in a square
for i=1:n_Emitters
    Emitters(i).HTM =Trans3(temp_x(i),temp_y(i),0)* Em_Base_HTM*...
        RotY3(pi/8*0);
end

%% Create the receivers

% Receivers:
Np = params.Np;                 % Number of parallels in the sensor
Nm = params.Nm;                 % Number of meridians in the sensor
n_Receivers = Np*Nm;    % Number of receivers
Ar = 1e-6;              % Active receiving area
Ts = 1;                 % Optical filter gain
n = 1;                  % Receiver's internal refractive index
Psi = params.Psi*pi/180;             % Hemi-Fov
R = 1;                  % Receiver's responsivity

% Create the receiver structure:
Receivers = newReceivers(n_Receivers,Ar, Ts, n, Psi, R);
% Receivers are organized in Parallel and Meridians arragement of photo
% detectors, with Nm Meridians and 3 Parallels, in a sphere with
% radius SR
SR = params.SR;
PDSensor = vlpCreateSensorParMer(Receivers, Np, Nm, SR, pi/8*0);

%% setup the values for computing received indication

% Quantities required for computation
%       Bw -      Bandwidth of receiver circuit
%       Z -       Vector with the transimpedance feedback resistors
%       s_i -     Vector with the operational amplifiers current PSD
%       s_v -     Vector with the operational amplifiers voltage PSD
%       Z_p -     Vector with the photo-diode equivalent impedances
%       Theta -   Thermodynamic temperature of feedback resistor, in Kelvin

Bw = 10e4;      % Bandwidth= 10kHz
Theta = 273+30;   % Feedback resistor at 30 degrees C

% Vector of ones for the receivers
nRec_v = ones(n_Receivers,1);

s_i = 1.3e-15*nRec_v;    % Current noise "plateau" at 1.3 fA/sqrt(Hz)
s_v = 4.8e-9*nRec_v;     % Voltage noise "plateau" at 4.8 nV/sqrt(Hz)
Z = 1e6*nRec_v;          % Feedback resistors = 1M
Z_p = 100e6*nRec_v;      % PD equivalent impedace = 100 MOhm


%% Main cycle

Wstep =0.5;
Lstep = 0.5;

xloc = 0:Wstep:W;
yloc = 0:Lstep:L;

Nrep = 1;

tempSensor = PDSensor;

for ix = 1:numel(xloc)
    for iy = 1:numel(yloc)
        % Move the sensor
        % Apply the transformation to every HTM in the sensors
        for i = 1:numel(tempSensor)
            tempSensor(i).HTM = Trans3(xloc(ix),yloc(iy),0)*PDSensor(i).HTM;
        end
        
        
        %  Compute received indication (mean and noise / variance)
        [ Y, nu ] = vlpRecIndication( Emitters, tempSensor, Bw, Z, s_i, s_v, Z_p, Theta );
        Nu = repmat(nu,1,n_Emitters);
        
        for counter = 1:Nrep
            s = sqrt(Nu).*randn(size(Y));
            Ynoise = Y ;%+ s;
            
            %             stemp(counter,:,:) = s;
            
            % Get a matrix with all HTMs, side-by-side
            x = [tempSensor.HTM];
            E = x(1:3,3:4:end);
            
            % Mvec is a matrix with the vectors pointing to the light sources
            Mvec = E*Ynoise;
            % Normalize Mvec
            Mvec = Mvec./repmat(sqrt(sum(Mvec.^2)),3,1);
            
            % The angle with the vertical is given by acos(kz*Mvec), where
            % kz = [0 0 1] (a vector pointing up). The internal product is simply
            % the third line of Mvec. The norm of both vectors is 1 (Mvec has
            % been normalized), so the expression can be simplified.
            vangles = acos(Mvec(3,:));
            
            % Compute the distances to light sources in the xy plane
            radii = H*tan(vangles);
            
            % recP holds the total power received from each emitter
            recP = sqrt(sum(Ynoise.*Ynoise));
            
            % Criteria for accepting the location data
            accepted=zeros(1,n_Emitters);
            mrecP = mean(recP);
            count = 0;
            while(sum(accepted) <3 && count< 300)
                count = count+1;
                accepted = recP > mrecP;
                mrecP = 0.98*mrecP;
            end
            
            if(count >= 300)
                location= [inf; inf];
                res(ix,iy).loc = location;
                res(ix,iy).pos = [xloc(ix);yloc(iy)];
                res(ix,iy).dist= norm(location-[xloc(ix);yloc(iy)]);
                res(ix,iy).radii = 0;
                res(ix,iy).accepted = accepted;
                res(ix,iy).trueradii= 0;
                res(ix,iy).broken = 1;
                broken=1;
                break;
            end
            
            % Get the emitters position
            temp = [Emitters.HTM];
            posEm = temp(1:2,4:4:end);
            %Consider only accepted positions
            posEm_ac = posEm(:,accepted);
            
            % Consider only accepted radii
            radii_ac = radii(accepted);
            
            Xx = posEm_ac(1,:);
            Yy = posEm_ac(2,:);
            
            A=[ Xx(2:end)'-Xx(1) Yy(2:end)'-Yy(1)];
            B=0.5*((radii_ac(1)^2-radii_ac(2:end)'.^2) + (Xx(2:end)'.^2+Yy(2:end)'.^2) - (Xx(1)^2 + Yy(1)^2));
            
            location = (A'*A)^(-1)*A'*B;
            
            % Compute the true value for radii
            delta = posEm - repmat([xloc(ix);yloc(iy)],1,n_Emitters);
            trueradii = sqrt(sum(delta.^2));
            
            %Send data to structure to save
            res(ix,iy).loc = location;
            res(ix,iy).pos = [xloc(ix);yloc(iy)];
            res(ix,iy).dist= norm(location-[xloc(ix);yloc(iy)]);
            res(ix,iy).radii = radii;
            res(ix,iy).accepted = accepted;
            res(ix,iy).trueradii= trueradii;
        end
    end
end


%calculate the percentage of values under 10 cm
max_H = 0.1;
mtemp=reshape([res.dist],ix,iy);
mtempb=double(mtemp < max_H);
mtemp=mtemp.*(1*mtempb-1e-2);
underPercentage=sum(sum(mtempb))/numel(mtemp)

