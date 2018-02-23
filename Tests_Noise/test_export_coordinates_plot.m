%%
experiment_name = 'WACOWCmulti2';
testWACOWCmulti(experiment_name);
clearvars -except experiment_name;
run extract_testWACOWCmulti.m

%%
experiment_number=  2;%1 to 4
test_number = 15 ; % 1 to 65 -> information about the test in the params structure


location_xy_plane = export_coordinates( experiment_number, test_number, 5);

%%


load(['results/data_WACOWCmulti' num2str(experiment_number) '.mat']);

temp_data = data(test_number);
% Params Structure
params = temp_data.params;

x_index= 15;
y_index = 15;

locations=location_xy_plane(x_index,y_index).locations;





% test_data.export(x_index, y_index).locations;
% 
% 
% 
% % for easier handling you can do
% 
% locations = test_data.export(x_index,y_index).locations;

%this way you direct access to the array of locations



%for visualization purposes

plot(locations(1,:),locations(2,:),'.')

hold on

axis([0 4 0 4])

plot(params.Wstep*(x_index-1),params.Lstep*(y_index-1), 'or')

