function PlotEmittersAndIntensity( Emitters, xvals, yvals, vals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Plot the light sources HTM reference frames
for e = Emitters
    plot3Drefaxis(e.HTM)
    hold on
    axis equal
end

h = surf(xvals,yvals,vals);
h.MeshStyle = 'none';
% view(2);
% figure(gcf);
axis equal

end

