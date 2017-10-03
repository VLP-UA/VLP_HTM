function [ Emitters ] = newEmitters( varargin )
%NEWEMITTERS Creates new Emitters
%   Emitters = newEmitters(n_Emitters, Pt, m);
%   Emitters = newEmitters(Pt, m);
%   Emitters = newEmitters(n_Emitters);
%   Emitters = newEmitters();
%   
%   n_Emitters  : number of emitters to create
%   Pt          : transmitted power
%   m           : Lambertian mode number
%
%   Creates an emitter of an array of n_Emitters emitters. If n_Emitters is
%   omitted, Emitters has 1 single element. If Pt and m are omitted, Pt
%   defaults to 0 and m defaults to 1. 
%
%   The HTM of the emitters are placed at the origin.


narginchk(0,3);

if (nargin == 3)
    n_Emitters = varargin{1};
    Pt = varargin{2};
    m = varargin{3};
elseif (nargin == 2)
    n_Emitters = 0;
    Pt = varargin{1};
    m = varargin{2};
elseif (nargin == 1)
    n_Emitters = varargin{1};
    Pt = 0;
    m = 1;
else
    n_Emitters = 0;
    Pt = [];
    m = [];
end

% Create the emitter structure:
Emitter_t = struct('HTM',{},'Pt',{},'m',{});

if(numel(Pt) == 1)
    % Pt and m are defined
    Emitter_t(1).Pt = Pt;
    Emitter_t(1).m  = m;
    Emitter_t(1).HTM = eye(4);
end


if (n_Emitters == 0)
    % Emitters is just a single struct
    Emitters = Emitter_t;
    return;
else
    % Replicate to create the receivers array:
    Emitters = repmat(Emitter_t,1,n_Emitters);
end

% TODO: Pt and m are [n_Emitters x 1] arrays with the individual values of
% Pt and m for every emitter

end

