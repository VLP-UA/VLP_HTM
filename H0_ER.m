function [ H ] = H0_ER( HTM_E, HTM_R, m, Ar, varargin )
% H0_ER( HTM_E, HTM_R, m, Ar, @g, @T ) 
%   OWC gain 
%   Computes the DC gain of a optical channel from emitter to receiver
%   following the models defined by Gfeller and Bapst [1] and Barry et al
%   [2]. 
%
%   The position and orientation of emitter and receiver are defined
%   by Homogeneous Transformation Matrices HTM_E and HTM_R. The axis of
%   emitter and receiver are aligned with the respective z axis.
%
%   m is the Lambertian order of the emitter, Ar the sensor area at the
%   receiver, g and T are, respectively, the gain and filter at the
%   receiver.
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

% Compute distance from emitter to receiver
d = norm(P_e - P_r);

% Compute u, the versor in the direction P_e to P_r: 
u = (P_r - P_e)/d;

if nargin > 4 
    % g and T functions present in the argument list
    % incidence angle must be computed
    Theta = acos(dot(n_r,-u));
    gf = varargin{1};
    g = gf(Theta);
    if nargin == 6
        Tf = varargin{2};
        T = Tf(Theta);
    else
        T = 1;
    end
else
    g = 1; 
    T = 1; 
end

H = (m+1)/(2*pi)*Ar*T*g*dot(n_e,u)*dot(n_r,-u)/d^2;

end

