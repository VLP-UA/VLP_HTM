
clear all;

% add the path to the projective geometry functions
addpath('../../ProjGeom');
addpath('../');



%% simulation parameters
%emitter
params.n_Emitters = 4; %square mode
params.m =3;
params.Pb = 1;
params.Ps = 0.25;

params.Re = 2.5;%obsolete

%room
params.W = 7.5;%not to be changed
params.L = 7.5;%not to be changed
params.H = 2.5;%not to be changed

%receiver
params.Np = 18;
params.Nm = 23;
params.Psi = 20;
params.SR = 0.25;




SR1=[0.25 0.2 0.15];
m1=[1,2,3,4,5,6];
Psi1=[15 18 22.5 30];
Np1=[6:20];
Nm1= [11:25];

for a=1:numel(SR1)
    
    for b=1:numel(m1)
        
        for c=1:numel(Psi1)
            
            for d=1:numel(Np1)
                
                 for e=1:numel(Nm1)
                     
                     params.SR = SR1(a);
                     params.m =  m1(b);
                     params.Psi = Psi1(c);
                     params.Np = Np1(d);
                     params.Nm = Nm1(e);
                     
                     run run_mult_parameters_daugther.m
                     clearvars -except a b c d e m1 Np1 SR1 Nm1 Psi1 params
                 end
            end
        end
    end
end




