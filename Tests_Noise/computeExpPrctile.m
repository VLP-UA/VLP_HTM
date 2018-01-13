function [ prctil_val, Psi_v, Np_v ] = computeExpPrctile( experiment_name , p, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Get the list of all mat files with test data
list = dir([experiment_name '_1*.mat']);

if numel(list) == 0
  error('No files found!\n');
end

% Get values in first file to initialize storage
load (list(1).name );

Psi_v = params.Psi_v(:);
Np_v = params.NmNp(:,2);

Psi_number = numel(Psi_v);
Np_number = numel(Np_v);

% Create the array
prctil_val = zeros(Psi_number,Np_number);
prctil_val_restrict = zeros(Psi_number,Np_number);

% Iterate over all test files
for i = 1:numel(list)
  matfilename = list(i).name;
  
  load(matfilename);
  
  % Compute Psi_min and Psi_r (ratio of Psi to Psi_min)
  Psi_min = 0.5*acos(cos(pi./(2*params.Np)).^2);
  Psi_r = round(params.Psi/Psi_min,2);
  
  % Get the indexes for Psi and Np values
  Psi_index = find(params.Psi_v == Psi_r);
  Np_index = find(Np_v==params.Np);
  
  % Store percentile value in array
  prctil_val(Psi_index, Np_index) = prctile(results.locerrorrms(:), p);
  
  % Restricted ared indexes:
  [npW npL] = size(results.locerrorrms);
  rl = round([1/6 5/6]' * [npW  npL]) ;
  
  prctil_val_restrict(Psi_index, Np_index) = ...
    prctile(results.locerrorrms(rl(1,1):rl(1,2),rl(2,1):rl(2,2)), p);
  
end

% Handle legacy error on variable name
if isfield(params,'N_L')
  params.n_L = params.N_L;
  params.n_W = params.N_W;
end

plotGraphics = 1;
map_threshold = 1;

if(plotGraphics)
  figure
  contourf(Psi_v,Np_v,prctil_val');
  colorbar
  title(['rms error, ' num2str(p) '-percentile (' ...
    num2str(params.n_W) 'x' num2str(params.n_L) ') light sources' ]);
  xlabel('\Psi_r')
  ylabel('N_p')
  caxis([0 map_threshold])
  map = colormap;
  map(size(map,1),:) = [1 1 1];
  colormap(map);

  figure
  contourf(Psi_v,Np_v,prctil_val_restrict');
  colorbar
  title(['rms error, ' num2str(p) '-percentile (' ...
    num2str(params.n_W) 'x' num2str(params.n_L) ') light sources, restricted area' ]);
  xlabel('\Psi_r')
  ylabel('N_p')
  caxis([0 map_threshold])
  map = colormap;
  map(size(map,1),:) = [1 1 1];
  colormap(map);

end

% Report results
disp(sprintf('%s: ',experiment_name))
disp(sprintf('Full area:  rms: %f min: %f',rms(prctil_val(:)),min(prctil_val(:))))
disp(sprintf('Restricted: rms: %f min: %f\n',rms(prctil_val_restrict(:)),min(prctil_val_restrict(:))))


end

