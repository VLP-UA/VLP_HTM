function [ h ] = drawCircleAtH( Cx, Cy, H, R )
%DRAWDIRCLEATH Draw a circle at height H, centered in Cx,Cy with radius R
%   Detailed explanation goes here

theta = linspace(0,2*pi);

X = Cx + R*cos(theta);
Y = Cy + R*sin(theta);
Z = H*ones(size(X));

h = plot3(X,Y,Z,'b--');

end

