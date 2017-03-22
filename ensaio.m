clear t 

delta = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001]

for ind=1:numel(delta)
    
    clearvars -except t delta dd ind
    
    % Define the delta increment
    dx = delta(ind)
    dy = delta(ind);
    
    DefineSpace
    
    
    % 1st method
    
    tic
    experiment3
    t(1,ind)=toc
    
    %figure
    %PlotEmittersAndIntensity(Emitters,xvals,yvals,rec_power);
    
    % 2nd method
    
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

figure
loglog((Lx./delta).^2,t(1,:)./t(2,:),'*-')
xlabel('# points')
ylabel('Ratio of computation time (Meth 1/Meth 2)')
title('Computation time for method 1 and 2')

s = [ 'ComparisonResult_' datestr(now,'yyyymmdd_hhMMss')];
save(s, 'delta', 't');
