clear t 

delta = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001]

for ind=1:numel(delta)
    
    clearvars -except t delta dd ind
    
    % Define the delta increment
    dx = delta(ind)
    dy = delta(ind);
    
    DefineSpace
   
    CreateReceiversWithHTM
    
    tic
    experiment_2_1
    t(2,ind)=toc
    
    %figure
    %PlotEmittersAndIntensity(Emitters,xvals,yvals,vals);
    
end


figure
plot((Lx./delta).^2,t,'*-')
xlabel('# points')
ylabel('Computation time (s)')
title('Computation time for method 1 and 2')
legend('Method 1','Method 2','Location','NorthWest')

figure
loglog((Lx./delta).^2,t,'*-')
xlabel('# points')
ylabel('Computation time (s)')
title('Computation time for method 1 and 2')
legend('Method 1','Method 2','Location','NorthWest')

s = [ 'ComparisonResult_' datestr(now,'yyyymmdd_hhMMss')];
save(s, 'delta', 't');
