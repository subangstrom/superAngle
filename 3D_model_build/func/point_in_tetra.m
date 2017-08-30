function [ chk ] = point_in_tetra( point, model_mesh_nei, tetra_num )
%Weizong Xu, August, 2017

tetra_P1=model_mesh_nei.p(:,model_mesh_nei.t_sort(1,tetra_num))';
tetra_P2=model_mesh_nei.p(:,model_mesh_nei.t_sort(2,tetra_num))';
tetra_P3=model_mesh_nei.p(:,model_mesh_nei.t_sort(3,tetra_num))';
tetra_P4=model_mesh_nei.p(:,model_mesh_nei.t_sort(4,tetra_num))';
tot_V = Vcal_tetrahedron( point,tetra_P1,tetra_P2,tetra_P3 ) + ...
        Vcal_tetrahedron( point,tetra_P1,tetra_P2,tetra_P4 ) +...
        Vcal_tetrahedron( point,tetra_P1,tetra_P3,tetra_P4 ) + ...
        Vcal_tetrahedron( point,tetra_P2,tetra_P3,tetra_P4 );
tetra_V=model_mesh_nei.tetra_volume{tetra_num,1};

if ((tot_V-tetra_V)<=1e-10)
    chk=1; %it is in tetrahedral
else
    chk=0;
end

end

