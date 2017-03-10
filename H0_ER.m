function [ H ] = H0_ER( HTM_E, HTM_R, m, Ar, varargin )
% H0_ER( HTM_E, HTM_R, m, Ar ) 
% H0_ER( HTM_E, HTM_R, m, Ar, FOV ) 
% H0_ER( HTM_E, HTM_R, m, Ar, @g ) 
% H0_ER( HTM_E, HTM_R, m, Ar, @g, @T ) 
%
%   OWC channel gain 
%
%   Computes the DC gain of a optical channel from emitter to receiver
%   following the models defined by Gfeller and Bapst [1] and Barry et al
%   [2]. 
%
%   The position and orientation of emitter and receiver are defined
%   by Homogeneous Transformation Matrices HTM_E and HTM_R. The axis of
%   emitter and receiver are aligned with the respective z axis.
%
%   m is the Lambertian order of the emitter, Ar the sensor area at the
%   receiver; FOV is the semi-angle for the Field of View; g and T are,
%   respectively, the gain and filter at the receiver.
%
%   If g and T are ommitted, they default to 1. 
%
%   HTM_E : HTM with position and orientation of emitter;
%   HTM_R : HTM with position and orientation of receiver;
%   m : Lambertian mode number of emitter; 
%   Ar : area of receiver detection surface;
%   g : gain of optical concentrator;
%   T : optical filter transfer function
%
%   [1] F. R. Gfeller and U. Bapst, ‘Wireless in-house data communication 
%   via diffuse infrared radiation’, Proceedings of the IEEE, vol. 67, no. 
%   11, pp. 1474–1486, 1979.
%   [2] J. R. Barry, J. M. Kahn, W. J. Krause, E. A. Lee, and D. G. 
%   Messerschmitt, ‘Simulation of multipath impulse response for indoor 
%   wireless optical channels’, IEEE Journal on Selected Areas in 
%   Communications, vol. 11, no. 3, pp. 367–379, Apr. 1993.

% TODO : acrescentar FOV

% Get emiiter and receiver position
P_e = HTM_E(1:3,4);
P_r = HTM_R(1:3,4);

% Get emitter and receiver orientation. 
n_e = HTM_E(1:3,3);
n_r = HTM_R(1:3,3);

% Receiver and emitter must be facing each other
if dot(n_e,n_r) > 0
    H = 0;
    return
end

% Compute distance from emitter to receiver
d = norm(P_e - P_r);

% Compute u, the versor in the direction P_e to P_r: 
u = (P_r - P_e)/d;

% Check 5th argument
if nargin > 4 
    % incidence angle must be computed
    Theta = acos(dot(n_r,-u));

    gf = varargin{1};
    % Check if 5th argument is a function handle or a numerical value
    if isa(gf,'function_handle') 
        % Gain is defined by a function
        g = gf(Theta);
    else
        % A FOV was defined
        g = rectangularPulse(-gf,gf,Theta);
    end

    % Check 6th argument
    if nargin == 6
        Tf = varargin{2};
        if isa(Tf,'function_handle')
            % Filter is defined by a function
            T = Tf(Theta);
        else
            % The filter is defined by a constant
            T = Tf;
        end
    else
        T = 1;
    end
else
    g = 1; 
    T = 1; 
end

H = (m+1)/(2*pi)*Ar*T*g*(dot(n_e,u)^m)*dot(n_r,-u)/d^2;

if H<0
    H=0;
end

end

