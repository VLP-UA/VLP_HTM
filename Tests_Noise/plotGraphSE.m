function [hh] = plotGraphSE( resultsdir, resultfn, m, Np, Nm, Psi, Nrep, varargin)
%PLOTGRAPHSE Plots the graphs from single emitter using data in file
%
% plotGraphSE( resultsdir, resultfn, m, Np, Nm, Psi, Nrep)
% plotGraphSE( resultsdir, resultfn, m, Np, Nm, Psi, Nrep,t)
%
% Plots the results generated by test_single_em

[res params] = shapeResults(resultsdir, resultfn,m,Np,Nm,Psi,Nrep);

% Position in x and y
xloc = 0:params.Wstep:params.W;
yloc = 0:params.Lstep:params.L;

baseTitle = [ '\{m, N_p, N_m, \Psi \}=\{' num2str(m) ',' ...
  num2str(Np) ',' num2str(Nm) ',' ...
  num2str(round(180/pi*Psi)) '\}'];

% By default, plot all graphs
plotmax = 1;    % Plot max
plotavg = 1;    % Plot mean
plotstd = 1;    % Plot std dev

if nargin >= 8
  t = varargin{1};
  if isstr(t)
    switch t
      case 'max'
        plotavg = 0;
        plotstd = 0;
      case 'avg'
        plotmax = 0;
        plotstd = 0;
      case 'std'
        plotmax = 0;
        plotavg = 0;
      otherwise
        warning('t must be one of [max|avg|std]');
    end
  else
    warning('t must be a string!');
  end
end

% Vector for storing the graphics handles
hh = zeros(6,1);

if plotmax
  % Location error,max
  figure
  hh(1) = surf(xloc,yloc,(res.locerrormax)');
  t = title(['Location, max error | ' baseTitle]);
  set(t,'FontWeight','Normal');
  
  % Radius error,max
  figure
  hh(2) = surf(xloc,yloc,(res.raderrormax)');
  t = title(['Radius, max error | ' baseTitle]);
  set(t,'FontWeight','Normal');
end

if plotavg
  % Location error,mean
  figure
  hh(3) = surf(xloc,yloc,(res.locerroravg)');
  t = title(['Location, mean error | ' baseTitle]);
  set(t,'FontWeight','Normal');
  
  % Radius error,mean
  figure
  hh(4) = surf(xloc,yloc,(res.raderroravg)');
  t = title(['Radius, mean error | ' baseTitle]);
  set(t,'FontWeight','Normal');
end

if plotstd
  % Location error, stddev
  figure
  hh(5) = surf(xloc,yloc,(res.locerrorstd)');
  t = title(['Location , st. dev.| ' baseTitle]);
  set(t,'FontWeight','Normal');
  
  % Radius error,stdev
  figure
  hh(6) = surf(xloc,yloc,(res.raderrorstd)');
  t = title(['Radius , st. dev.| ' baseTitle]);
  set(t,'FontWeight','Normal');
end

end

