function plotExpGraphSE1(experiment_name)

% Load the experiment parameters
% run(experiment_name);

%% Files and folders

baseFileName = [ experiment_name '_data.mat' ];
resultsdir = '';


if exist(baseFileName) ~= 2
  % if the file does not exist, try results/ subdir
  % Folter for storing results
  resultsdir = 'results/';
end


% Experimente results filename
expFilename = [ resultsdir baseFileName ];

load (expFilename);

Ns_v = prod(params.NmNp,2);
Psi_v = params.Psi_v(:);

% testresultstable = struct2table(testresults);

m_v = params.m_v;

plotBaseFn = [resultsdir experiment_name];
%% 
% Test: Psi
% Filter: m
% Implicit: Np,Nm
for m = m_v
  
  % filter the table on m
  ifx = testresultstable.m==m;
  filtered = testresultstable(ifx,:);
  
  % Rearrange Psi values
  % Psi values are a function of Np. Replace actual values by the ration
  % Psi/Psi_min
  % 1. Get the values for corresponding Phi_min
  Psi_min = 0.5*acos(cos(pi./(2*filtered.Np)).^2);
  % 2. Compute the ratio Psi/Psi_min
  Psi_ratio = round(filtered.Psi ./ Psi_min,2);
  % 3. Update table
  filtered.Psi = Psi_ratio;
  
  
  % Sort by Implicit and then by Test
  [x ix] = sort(filtered.Np);
  filtered = filtered(ix,:);
  [x ix] = sort(filtered.Psi);
  filtered = filtered(ix,:);
  
  % Get the indexes for the plots (unique values for each var)
  valstest = unique(filtered.Psi);
  valsimpl = unique(filtered.Np);
  
  nvalstest = numel(valstest);
  nvalsimpl = numel(valsimpl);
  
  % Plot location error
  % Create the array:
  plotdata = reshape(filtered.locerrrms, nvalsimpl,nvalstest);
  
  h = figure;
  plot(Psi_v,plotdata','*-')
  ylim([0 1.1])
  
  title(['location, rms error, m=' num2str(m) ', VRT=' num2str(params.validReadingThreshold) ])
  ylabel('rms error/m')
  xlabel('\Psi/\Psi_{min}')
  legend([ repmat('N_s=',size(Ns_v)) num2str(Ns_v) ]);
  
  filename = [ plotBaseFn '_LocRMSErr_Psi_m_' num2str(m) '.png' ];
  print(filename,'-dpng');


  if exist('filtered.raderrrms')
    % Plot radius error
    % Create the array:
    plotdata = reshape(filtered.raderrrms, nvalsimpl,nvalstest);
    
    
    h = figure;
    plot(Psi_v,plotdata','*-')
    ylim([0 Inf])
    title(['radius, rms error, m=' num2str(m) ', VRT=' num2str(params.validReadingThreshold) ])
    ylabel('rms error/m')
    xlabel('\Psi/\Psi_{min}')
    legend([ repmat('N_s=',size(params.NmNp,1),1) num2str(prod(params.NmNp,2)) ]);
    
    filename = [ plotBaseFn '_RadRMSErr_Psi_m_' num2str(m) '.png' ]
    print(filename,'-dpng')
    
  end
  
end

%%
% Test: Np,Nm
% Filter: m
% Implicit: Psi
for m = m_v
  
  % filter the table on m
  ifx = testresultstable.m==m;
  filtered = testresultstable(ifx,:);
  
  % Rearrange Psi values
  % Psi values are a function of Np. Replace actual values by the ratio
  % Psi/Psi_min
  % 1. Get the values for corresponding Phi_min
  Psi_min = 0.5*acos(cos(pi./(2*filtered.Np)).^2);
  % 2. Compute the ratio Psi/Psi_min
  Psi_ratio = round(filtered.Psi ./ Psi_min,2);
  % 3. Update table
  filtered.Psi = Psi_ratio;
  
  shortlist = Psi_ratio==1 | Psi_ratio==2 | Psi_ratio==3 | Psi_ratio==4 | Psi_ratio==5;
  
  filtered = filtered(shortlist,:);
  
  % Sort by Implicit and then by Test
  [x ix] = sort(filtered.Psi);
  filtered = filtered(ix,:);
  [x ix] = sort(filtered.Np);
  filtered = filtered(ix,:);
  
  % Get the indexes for the plots (unique values for each var)
  valstest = unique(filtered.Np);
  valsimpl = unique(filtered.Psi);
  
  nvalstest = numel(valstest);
  nvalsimpl = numel(valsimpl);
  
  % Plot location error
  % Create the array:
  plotdata = reshape(filtered.locerrrms, nvalsimpl,nvalstest);
  
  h = figure;
  plot(params.NmNp(:,2),plotdata','*-')
  ylim([0 1.1])
  
  % Handle legacy error on variable name
  if isfield(params,'N_L')
    params.n_L = params.N_L;
    params.n_W = params.N_W;
  end

  title(['Positioning r.m.s. error (' ...
    num2str(params.n_W) 'x' num2str(params.n_L) ') light sources' ]);
  ylabel('rms error/m')
  xlabel('# of parallels')
  legend([ repmat('\Psi_r=',size(valsimpl)) num2str(valsimpl) ]);
  grid on
  
  filename = [ plotBaseFn '_LocRMSErr_Ns_m_' num2str(m) '.png' ];
  print(filename,'-dpng')
  % createfigureWACOWCpos(params.NmNp(:,2),plotdata');

  if exist('filtered.raderrrms')
    % Plot radius error
    % Create the array:
    plotdata = reshape(filtered.raderrrms, nvalsimpl,nvalstest);
    
    
    h = figure;
    plot(Psi_v,plotdata','*-')
    ylim([0 Inf])
    title(['radius, rms error, m=' num2str(m) ', VRT=' num2str(params.validReadingThreshold) ])
    ylabel('rms error/m')
    xlabel('# of photo-detectors')
    legend([ repmat('\Psi_r=',size(Psi_v)) num2str(Psi_v) ]);
    
    filename = [ plotBaseFn '_RadRMSErr_Ns_m_' num2str(m) '.png' ]
    print(filename,'-dpng')
    
    
  end
end

end
