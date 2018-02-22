
function [location_xy_plane] = export_coordinates( experiment_number, test_number, Ng)
%% export_coordinates - export the data
% experiment_number - number of the experiment to run -> 1 : 4
% test_number - number of the test to run -> 1:65
% Ng - number of elements per group -> 3:params.nEmitters
%
% %


%% add the path to the projective geometry functions
addpath('../../ProjGeom');
addpath('../');

% clear all;
% clc;


%% export from raw data to xy coordinates

% experiment_number=  2;%1 to 4
% test_number = 50 ; % 1 to 65 -> information about the test in the params structure

load(['results/data_WACOWCmulti' num2str(experiment_number) '.mat'])

temp_data = data(test_number);
% Params Structure
params = temp_data.params;

% Grouping parameters
Ng = 3; % number of emitters in the group
use_average=0; % set to one to use the average of params.Nrep repetitions

%% generate emitters from params
%%Create the light emitters
Emitters = newEmitters(params.n_Emitters,params.Pb,params.Ps,params.m);
% Emitter location in xy plane
locEm = [ 0 0 params.H ];
Em_Base_HTM = Trans3( locEm' )*RotX3(pi);      % Base HTM at (0,0), on the ceiling.
delta_W = params.W/(params.n_W+1);
delta_L = params.L/(params.n_L+1);
for i=1:params.n_W
    for j = 1:params.n_L
        Emitters(j+(i-1)*params.n_L).HTM = Trans3(i*delta_W,j*delta_L,0)*Em_Base_HTM;
    end
end



%% Iterate over the xy plane and aggregate the emmiters into groups


for xi=1: size(temp_data.export_radii,1)
    for yi=1: size(temp_data.export_radii,2)
        
        %extract variables from structure for easier handling
        export = temp_data.export_radii(xi,yi);
        
        radii_all = export.computed_radii;
        power = export.rec_power;
        
        if(use_average)
            radii_all_ave = export.average_radii;
            power_average = export.average_rec_power;
        end
        
        %generate the combination of N
        combi= nchoosek([1:params.n_Emitters],Ng);
        
        
        loc_xy=[];
        loc_ave_xy = [];
        power_xy =[];
        power_xy_ave=[];
        % Calculate an estimate for each group of photodiodes
        for i =1:size(combi,1)
            radii = radii_all(combi(i,:));
            
            location = trilateration_group_em(Emitters(combi(i,:)), radii);
            loc_xy = [loc_xy location];
            power_xy = [power_xy sum(power(combi(i,:)))];
            
            if(use_average)
                radii_ave = radii_all_ave(combi(i,:));
                location_average = trilateration_group_em(Emitters(combi(i,:)), radii_ave);
                loc_ave_xy = [loc_ave_xy location_average];
                power_xy_ave = [power_xy_ave sum(power_average(combi(i,:))) ];
            end
            
            
        end
        
        location_xy_plane(xi,yi).locations = loc_xy;
        location_xy_plane(xi,yi).power = power_xy;
        
        if(use_average)
            location_xy_plane(xi,yi).power_ave= power_xy_ave;
            location_xy_plane(xi,yi).location_ave = loc_ave_xy;
        end
    end
end
end

