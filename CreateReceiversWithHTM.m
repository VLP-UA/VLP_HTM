%% Receivers
% The receiving surface (the floor) is emulated by a series of receivers,
% with FOV=pi/2 and no gain.

% Create the receiver structure:
Receiver_t = struct('HTM',{},'Ar',{},'Ts',{},'n',{},'Psi',{},'Pr',{});

% Create one element with default values
Receiver_t(1).Ar = dx*dy;
Receiver_t(1).Ts = 1;
Receiver_t(1).n = 1;
Receiver_t(1).Psi = pi/2;
Receiver_t(1).Pr = 0;

% Replicate to create the receivers array:
Receivers = repmat(Receiver_t,1,nx*ny);

% First receiver is at (0,0,0), looking up (no change in the axis
% orientation)
HTM_R_base = Trans3(0,0,0);


for ix = 1:nx
    for iy = 1:ny
        % Compute receiver HTM
        ri = iy + nx*(ix-1);
        Receivers(ri).HTM = Trans3(xvals(ix),yvals(iy),0);
    end
end
