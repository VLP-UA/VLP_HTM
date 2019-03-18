clear
close all
clc

%%
load('..\Clustering\epsilon_minPts_matrix_62_exp_4.mat')


%%
load('dbscan_data.mat');
load('kmeans_data.mat');
load('meanshift_data.mat');

%%

% add the path to the projective geometry functions
addpath('../../ProjGeom');
addpath('../');
addpath('../Tests_Noise\')

addpath('../Clustering\')


%%
first_n_groups_array = 5:10:55;
%%

for first_n_groups = first_n_groups_array
    
    plotOn = 0;
    group_overlap = 1;
    
    Nm = params.Nm;
    Np = params.Np;
    
    Psi=params.Psi;
    
    run partial_sim.m
    save(['data_first_n_groups_' num2str(first_n_groups)], 'ground');
    clearvars -except params first_n_groups DBSCAN_data KMEANS_data MEANSHIFT_data
end


%%
clearvars -except first_n_groups_array

corner_x = 20:41;
corner_y = 20:41;

center_x = 15:25;
center_y = 15:25;

index = 1;

for first_n_groups = first_n_groups_array
    load(['data_first_n_groups_' num2str(first_n_groups)]);
    temp_ground = reshape([ground(:,:).trilat_error],size(ground,1), size(ground,2));
    trilat_mean_ground_corner(index) = mean(mean(temp_ground(corner_x, corner_y)));
    trilat_mean_ground_center(index) = mean(mean(temp_ground(center_x, center_y)));
    trilat_mean_ground(index) = mean(mean(temp_ground));
    
    %     temp_ground = reshape([ground(:,:).dbscan_error_real],size(ground,1), size(ground,2));
    %     dbscan_mean_ground_corner(index) = mean(mean(temp_ground(corner_x, corner_y)));
    %     dbscan_mean_ground_center(index) = mean(mean(temp_ground(center_x, center_y)));
    %     dbscan_mean_ground(index) = mean(mean(temp_ground));
    
    
    temp_ground = reshape([ground(:,:).kmeans_error_real],size(ground,1), size(ground,2));
    kmeans_mean_ground_corner(index) = mean(mean(temp_ground(corner_x, corner_y)));
    kmeans_mean_ground_center(index) = mean(mean(temp_ground(center_x, center_y)));
    kmeans_mean_ground(index) = mean(mean(temp_ground));
    
    temp_ground = reshape([ground(:,:).meanshift_error_real],size(ground,1), size(ground,2));
    meanshift_mean_ground_corner(index) = mean(mean(temp_ground(corner_x, corner_y)));
    meanshift_mean_ground_center(index) = mean(mean(temp_ground(center_x, center_y)));
    meanshift_mean_ground(index) = mean(mean(temp_ground));
    
    
    index = index + 1;
    
end


figure
hold on
plot(first_n_groups_array,trilat_mean_ground_corner,'-')
% plot(first_n_groups_array,dbscan_mean_ground_corner,'*-')
plot(first_n_groups_array,kmeans_mean_ground_corner,'+-')
plot(first_n_groups_array,meanshift_mean_ground_corner,'o-')

figure
hold on
plot(first_n_groups_array,trilat_mean_ground_center,'-')
% plot(first_n_groups_array,dbscan_mean_ground_center,'*-')
plot(first_n_groups_array,kmeans_mean_ground_center,'+-')
plot(first_n_groups_array,meanshift_mean_ground_center,'o-')

figure
hold on
plot(first_n_groups_array,trilat_mean_ground,'-')
% plot(first_n_groups_array,dbscan_mean_ground,'*-')
plot(first_n_groups_array,kmeans_mean_ground,'+-')
plot(first_n_groups_array,meanshift_mean_ground,'o-')













