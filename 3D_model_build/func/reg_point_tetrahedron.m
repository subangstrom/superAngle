function [ model_mesh_nei, p_tetra_num ] = reg_point_tetrahedron( point, model_mesh_nei, chk_reg )
%Weizong Xu, August, 2017
%chk_reg 1.write p_tetra_num in the model file, -1 no display no write
%other no write but display
num=0;
for ii=1:size(model_mesh_nei.t_sort,2)
    %ii
    tetra_P1=model_mesh_nei.p(:,model_mesh_nei.t_sort(1,ii))';
    tetra_P2=model_mesh_nei.p(:,model_mesh_nei.t_sort(2,ii))';
    tetra_P3=model_mesh_nei.p(:,model_mesh_nei.t_sort(3,ii))';
    tetra_P4=model_mesh_nei.p(:,model_mesh_nei.t_sort(4,ii))';
%     tot_V = Vcal_tetrahedron( point,tetra_P1,tetra_P2,tetra_P3 ) + ...
%         Vcal_tetrahedron( point,tetra_P1,tetra_P2,tetra_P4 ) +...
%         Vcal_tetrahedron( point,tetra_P1,tetra_P3,tetra_P4 ) + ...
%         Vcal_tetrahedron( point,tetra_P2,tetra_P3,tetra_P4 );

    P21 = tetra_P2-point;
    P31 = tetra_P3-point;
    P41 = tetra_P4-point;  
    volume= abs(P21(1)*(P31(2)*P41(3)-P31(3)*P41(2))-P21(2)*(P31(1)*P41(3)-P31(3)*P41(1))+P21(3)*(P31(1)*P41(2)-P31(2)*P41(1)));
    P21 = tetra_P1-point; 
    volume= volume+ abs(P21(1)*(P31(2)*P41(3)-P31(3)*P41(2))-P21(2)*(P31(1)*P41(3)-P31(3)*P41(1))+P21(3)*(P31(1)*P41(2)-P31(2)*P41(1)));
    P31 = tetra_P2-point; 
    volume= volume+ abs(P21(1)*(P31(2)*P41(3)-P31(3)*P41(2))-P21(2)*(P31(1)*P41(3)-P31(3)*P41(1))+P21(3)*(P31(1)*P41(2)-P31(2)*P41(1)));
    P41 = tetra_P3-point; 
    volume= volume+ abs(P21(1)*(P31(2)*P41(3)-P31(3)*P41(2))-P21(2)*(P31(1)*P41(3)-P31(3)*P41(1))+P21(3)*(P31(1)*P41(2)-P31(2)*P41(1)));
    tot_V=1/6*volume;
    
    tetra_V=model_mesh_nei.tetra_volume{ii,1};
    if ((tot_V-tetra_V)<=1e-10)
        p_tetra_num=ii;
        num=num+1;
        if (chk_reg~=-1)
            %disp(['Find tetrahedron #', num2str(ii)])
        end
        %look_p=ii;
        %look_model_tetra( model_mesh_nei, look_p, point' );
    end
end

if (num==1 && chk_reg==1)
    model_mesh_nei.start_point=point;
    model_mesh_nei.start_tetra_num=p_tetra_num;
end

if (num>1)
    disp('Error! Not a good start point, it share more than one tetrahedron')
    p_tetra_num=nan;
end

if (num==0)
    if (chk_reg~=-1)
        disp('Error! Model does not contain this point.')
    end
        p_tetra_num=nan;
end

end

