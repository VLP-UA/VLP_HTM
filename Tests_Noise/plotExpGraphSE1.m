function plotExpGraphSE1(experiment_data)

run(experiment_data);

%% Files and folders

% Folter for storing results
resultsdir = 'results/';

% Base filename
resultsBaseFn = experiment_data;

% Experimente results filename
expFilename = [ resultsdir resultsBaseFn '_data' ];


load (expFilename);

Ns_v = prod(params.NmNp,2);
Psi_v = params.Psi_v(:);

% testresultstable = struct2table(testresults);

plotBaseFn = [resultsdir resultsBaseFn];
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
  ylim([0 Inf])
  title(['location, rms error, m=' num2str(m) ', VRT=' num2str(params.validReadingThreshold) ])
  ylabel('rms error/m')
  xlabel('\Psi/\Psi_{min}')
  legend([ repmat('N_s=',size(Ns_v)) num2str(Ns_v) ]);
  
  filename = [ plotBaseFn '_LocRMSErr_Psi_m_' num2str(m) '.png' ]
  print(filename,'-dpng')

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
  plot(Ns_v,plotdata','*-')
  ylim([0 Inf])
  title(['location, rms error, m=' num2str(m) ', VRT=' num2str(params.validReadingThreshold) ])
  ylabel('rms error/m')
  xlabel('# of photo-detectors')
  legend([ repmat('\Psi_r=',size(Psi_v)) num2str(Psi_v) ]);
  
  filename = [ plotBaseFn '_LocRMSErr_Ns_m_' num2str(m) '.png' ]
  print(filename,'-dpng')
  
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
