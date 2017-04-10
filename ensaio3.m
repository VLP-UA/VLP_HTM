% ensaio3.m
%
% Verification of the time taken by each step of the computation.

clear t 

delta = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 ]

for ind=1:numel(delta)
    
    clearvars -except t delta dd ind
    
    % Define the delta increment
    dx = delta(ind)
    dy = delta(ind);
    
    DefineSpace
   
    CreateReceiversWithHTM
    
    tic
    A = Receiver2A(Receivers);
    t(1,ind)=toc

    B = Emitter2B(Emitters);
    t(2,ind)=toc

    F = EmitRec2F(Emitters, Receivers);
    t(3,ind)=toc

    Pr = A*F*B*[Emitters.Pt]';
    t(4,ind)=toc
    
%     vals = reshape(Pr,[nx ny]);
%     figure
%     PlotEmittersAndIntensity(Emitters,xvals,yvals,vals);
    
end

% to contains the time per each operation
to = [ t(1,:)
    diff(t)
    sum(to) ]

figure
plot((Lx./delta).^2,to,'*-')
xlabel('# points')
ylabel('Computation time (s)')
title('Computation time for method 1 and 2')
legend('A','B','F','Pr','Total','Location','NorthWest')

figure
loglog((Lx./delta).^2,to,'*-')
xlabel('# points')
ylabel('Computation time (s)')
title('Computation time for method 1 and 2')
legend('A','B','F','Pr','Total','Location','NorthWest')

% s = [ 'ComparisonResult_' datestr(now,'yyyymmdd_hhMMss')];
% save(s, 'delta', 't');
