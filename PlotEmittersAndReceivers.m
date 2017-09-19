function PlotEmittersAndIntensity( Emitters, Receivers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Plot the light sources HTM reference frames
for e = Emitters
    plot3Drefaxis(e.HTM)
    hold on
    axis equal
end

for e = Receivers
    plot3Drefaxis(e.HTM)
    hold on
    axis equal
end


end

