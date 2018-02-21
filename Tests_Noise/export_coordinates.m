%% export from raw data to xy coordinates

experiment_number=  2;%1 to 4

load(['results/data_WACOWCmulti' num2str(experiment_number) '.mat'])


temp_data = data(50);

params = temp_data.params;


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



for xi=1: size(temp_data.export_radii,1)
    for yi=1: size(temp_data.export_radii,2)
        
        export = temp_data.export_radii(xi,yi);

        radii_all = export.computed_radii;


        combi= nchoosek([1:numel(export.computed_radii)],9);

        loc_out=[];
        for i =1:size(combi,1)
            radii = radii_all(combi(i,:));

            % Get the emitters position
            temp = [Emitters.HTM];
            pos_Em = temp(1:2,4:4:end);
            %Consider only accepted positions
            pos_Em=pos_Em(:,combi(i,:));
            % pos_Em = posEm(:,accepted);

            Xx = pos_Em(1,:);
            Yy = pos_Em(2,:);

            A=[ Xx(2:end)'-Xx(1) Yy(2:end)'-Yy(1)];
            B=0.5*((radii(1)^2-radii(2:end)'.^2) + (Xx(2:end)'.^2+Yy(2:end)'.^2) - (Xx(1)^2 + Yy(1)^2));

            location = (A'*A)^(-1)*A'*B;
            
            loc_out = [loc_out location];
            scattered_data(xi,yi).locations = loc_out;
        end
    end
    
end


% plot(loc_out(1,:),loc_out(2,:),'*');







