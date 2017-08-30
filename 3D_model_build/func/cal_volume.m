function [ model_mesh_nei ] = cal_volume( model_mesh_nei)
%Calculate volume of tetrahedral mesh
%Weizong Xu, August, 2017

for ii=1:size(model_mesh_nei.t_sort,2)
    [ Volume ] = mesh_volume_tetra(ii, model_mesh_nei );
    model_mesh_nei.tetra_volume{ii,1}=Volume;
end

v_tot=0;
for mesh_num=1:size(model_mesh_nei.t_sort,2)
Volume = model_mesh_nei.tetra_volume{mesh_num,1};
v_tot=v_tot+Volume;
end
disp(['Total volume is ', num2str(v_tot)]);

end

