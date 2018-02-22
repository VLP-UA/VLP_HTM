%% Data Export - sequence


experiment_name = 'WACOWCmulti1';
run extract_testWACOWCmulti.m

experiment_name = 'WACOWCmulti2';
run extract_testWACOWCmulti.m

experiment_name = 'WACOWCmulti3';
run extract_testWACOWCmulti.m

experiment_name = 'WACOWCmulti4';
run extract_testWACOWCmulti.m


%%each element of the data_xy is a point on the xy plane with the radius
%%for the combinations of Ng emitters

data_xy = export_coordinates(2, 50, 3);

