function PlotHTMArray( HTMarray )
%PLOTHTMARRAY Plots the HTMs in a struct array
%   PLOTHTMARRAY(HTMarray) plots the HTM in the struct array HTMarray.
%
%   HTMarray is a struct array where the struct contains a field, called
%   HTM, a 4x4 array representing a Homogeneous Transformation Matrix.
%   
%   X and y axes are plotted in blue; z axis is plotted in red.
%
%   See also PLOT3DREFAXIS

%   pf@ua.pt
    

for e = HTMarray
    plot3Drefaxis(e.HTM)
    hold on
    axis equal
end

end

