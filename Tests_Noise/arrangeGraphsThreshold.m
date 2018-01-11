function [ output_args ] = arrangeGraphsThreshold( hh, t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% For all graphs in hh
%for f = hh(:)
for i = 1:numel(hh)
  f = hh(i);
  ax = ancestor(f,'axes');
  if( isempty(ax) == 0)
    colorbar(ax);
    caxis(ax,[0 t]);
    map = colormap;
    map(size(map,1),:) = [1 1 1];
    colormap(ax,map);
    view(ax,2);
  end
end


end

