function [ Volume ] = mesh_volume_tetra( mesh_num, model_mesh )
%Weizong Xu, August, 2017

nP1=model_mesh.t(1,mesh_num);
nP2=model_mesh.t(2,mesh_num);
nP3=model_mesh.t(3,mesh_num);
nP4=model_mesh.t(4,mesh_num);
P1=model_mesh.p(:,nP1)';
P2=model_mesh.p(:,nP2)';
P3=model_mesh.p(:,nP3)';
P4=model_mesh.p(:,nP4)';
[ Volume ] = Vcal_tetrahedron( P1,P2,P3,P4 );

end

