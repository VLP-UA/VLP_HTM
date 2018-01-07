% This script compiles all data into a struture in order to produce the
% desired graph output, there are 2 static paramters, m and Psi, being the
% Np and Nm the parameter under test
% Each run for a given pair of m and Psi will have 26 values for Np and 63
% for Nm, a total of 1638 simulations

%% Basic setup
clear all;
close all;

% add the path to the projective geometry functions
addpath('../../ProjGeom');
addpath('../');



%% simulation parameters
%emitter
params.n_Emitters = 4; %square mode
params.m =1;
params.Pb = 1;
params.Ps = 0.25;

params.Re = 2.5;%obsolete

%room
params.W = 7.5;%not to be changed
params.L = 7.5;%not to be changed
params.H = 2.5;%not to be changed

%receiver
params.Np = 28;
params.Nm = 70;
params.Psi = 20;
params.SR = 0.25;

%simulation resolution
Wstep = 0.5;

%% Parameter to be changed in the simulation
% m1 = [1,2,3,4,5,6];
% Psi1 = [5 7.5 10 12.5];
M1 = [1 2 3 4 5 6];
Psi1 = [5 7.5 10 12.5];

%from the comment above choose a value to set m and Psi
% params.m = ;
% params.Psi = ;

%% Compile data into data struture
index=1;


for iPsi = 1:numel(Psi1)
    for iM = 1:numel(M1)
        
        params.m=M1(iM);
        params.Psi=Psi1(iPsi);
        
        %build filename to load
        field = fieldnames(params);
        
        filename= [];
        for i = [9 10 11 2]
            filename = [ filename num2str(getfield(params,field{i})) '-'] ;
        end
        load(['./res/square_' filename num2str(Wstep) '.mat'], 'res',...
            'params','underPercentage','broken');
        
        %place data from file in struct for easier access
        data(index).res=res;
        data(index).params=params;
        data(index).underPercentage=underPercentage;
        data(index).broken=broken;
        index=index+1;
    end
end


%%

%get maximum mean and standard deviation from the collected data
for index = 1: numel(data)
    data(index).max = max([data(index).res.dist]);
    
    data(index).mean = mean([data(index).res.dist]);
    data(index).std = std([data(index).res.dist]);
    
    iPsi=find(Psi1==data(index).params.Psi);
    iM=find(M1==data(index).params.m);
    
    mMax(iM,iPsi) = data(index).max;
    mMean(iM,iPsi) = data(index).mean;
    mStd(iM,iPsi) = data(index).std;
end

%% plot max,mean and std vs m
for i = 1:numel(M1)
    figure(1)
    plot(Psi1,mMax(i,:),'*-')
    hold on
    figure(2)
    plot(Psi1,mMean(i,:),'*-')
    hold on
    figure(3)
    plot(Psi1,mStd(i,:),'*-')
    hold on
end
figure(1)
title({'Max Error vs Psi';...
    [num2str(params.Np) ' Parallels, ' num2str(params.Nm) ' Meridians']})
figure(2)
title({'Mean Error vs Psi';...
    [num2str(params.Np) ' Parallels, ' num2str(params.Nm) ' Meridians']})
figure(3)
title({'Standard Deviation Error vs Psi';...
    [num2str(params.Np) ' Parallels, ' num2str(params.Nm) ' Meridians']})

for i = 1:3
    figure(i)
    
    set(gca,'XTick',Psi1)
    legend('m=1','m=2', 'm=3', 'm=4', 'm=5', 'm=6','Location','SouthWest')
    ylabel('Error(meters)')
    xlabel('Psi')
end

%% plot max,mean,std vs m
for i = 1:numel(Psi1)
    figure(4)
    plot(M1,mMax(:,i),'-*')
    hold on
    figure(5)
    plot(M1,mMean(:,i),'-*')
    hold on
    figure(6)
    plot(M1,mStd(:,i),'-*')
    hold on
end
figure(4)
title({'Max Error vs Lambertian Index';...
    [num2str(params.Np) ' Parallels, ' num2str(params.Nm) ' Meridians']})
figure(5)
title({'Mean Error vs Lambertian Index';...
    [num2str(params.Np) ' Parallels, ' num2str(params.Nm) ' Meridians']})
figure(6)
title({'Standard Deviation Error vs Lambertian Index';...
    [num2str(params.Np) ' Parallels, ' num2str(params.Nm) ' Meridians']})

for i = 4:6
    figure(i)
    
    set(gca,'XTick',(M1(1):1:M1(end)))
    legend('Psi = 5','Psi = 7.5', 'Psi = 10', 'Psi = 12.5')
    ylabel('Error(meters)')
    xlabel('m - Lambertian Index')
end





