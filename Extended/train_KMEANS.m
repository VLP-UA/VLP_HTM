clc;
clear;
close all;

%% Load Data

addpath('../Clustering/');
load coordinate_data.mat

plotOn=0;

n_max_clusters = 15;


cluster_array = 1: n_max_clusters;
%% Define the point under test

x_index_array = 1:size(ground,1);
y_index_array = 1:size(ground,2);


epsilon = linspace(0.001,0.3,25);%0.001:0.001:0.75;%linspace(0.001,1,50);%[ 0.2 0.1 0.05 0.025 0.01 0.005 0.0025]; % min distance between points in meters
MinPts =3;%[3:50];%[20 15 12 10 7 5 3 2]; %min size of the cluster

%% load coordinates from the data file
for x_index = x_index_array
    for y_index = y_index_array
        display(['X: '  num2str(x_index) ' Y:' num2str(y_index)]);
        coordinates=ground(x_index,y_index).coordinates;
        
        real_pos = ground(x_index,y_index).xy;
        estimated_pos = ground(x_index,y_index).estimated_location;
        
        %% Run DBSCAN Clustering Algorithm
        
        for n_clusters = cluster_array
            opts = statset('Display','off');
            [cidx, ctrs] = kmeans(coordinates, n_clusters, 'Distance','city', ...
                'Replicates',5, 'Options',opts);
            
            if plotOn
                figure;
                
                hold on;
                plot(ctrs(:,1),ctrs(:,2),'gx', 'MarkerSize', 12);
                legend toggle
                axis([0 4 0 4]);
                
                plot(real_pos(1), real_pos(2),'k*')
                plot(estimated_pos(1), estimated_pos(2),'m*')
            end
            real_error(n_clusters) = inf;
            for i = 1:n_clusters
                if plotOn
                    plot(coordinates(cidx==i,1),coordinates(cidx==i,2),'o');
                    hold on;
                end
                
                % find min error of all clusters
                temp_error = norm(ctrs(i)-real_pos);
                if(temp_error <= real_error(n_clusters))
                    real_error(n_clusters)=temp_error;
                end  
                
            end 
        end
        KMEANS_data(x_index, y_index).real_error = real_error;
        KMEANS_data(x_index, y_index).n_clusters = ...
            cluster_array(min(find(real_error == min(real_error))));

    end
end
min(min(reshape([KMEANS_data(:).real_error],41,41)))
mean(mean(reshape([KMEANS_data(:).real_error],41,41)))
max(max(reshape([KMEANS_data(:).real_error],41,41)))

%%

figure
surf(1:41,1:41,reshape([KMEANS_data(:).n_clusters],41,41))

figure
surf(1:41,1:41,reshape([KMEANS_data(:).real_error],41,41))


save('dbscan_data', 'DBSCAN_data')



%%
% d_plot = cluster_data(x_index,y_index);
% % for epsilon_index = epsilon
% %     for MinPts_index = MinPts;
% smallest_error_matrix = reshape([d_plot.test.smallest_error],numel(epsilon),numel(MinPts));
% contour3(MinPts,epsilon,(smallest_error_matrix),50)
% % shading interp
%
% %     end
% % end



