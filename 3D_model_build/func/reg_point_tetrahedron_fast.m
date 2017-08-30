function [ p_tetra_num ] = reg_point_tetrahedron_fast( list, point, p, t_sort, tetra_volume )
%Weizong Xu, August, 2017
%chk_reg 1.write p_tetra_num in the model file, -1 no display no write
%other no write but display

p_tetra_num=0; %set default 0, not find

for i=1:length(list)
    ii=list(i);
    tetra_P1=p(:,t_sort(1,ii))';
    tetra_P2=p(:,t_sort(2,ii))';
    tetra_P3=p(:,t_sort(3,ii))';
    tetra_P4=p(:,t_sort(4,ii))';
    tot_V = Vcal_tetrahedron( point,tetra_P1,tetra_P2,tetra_P3 ) + ...
        Vcal_tetrahedron( point,tetra_P1,tetra_P2,tetra_P4 ) +...
        Vcal_tetrahedron( point,tetra_P1,tetra_P3,tetra_P4 ) + ...
        Vcal_tetrahedron( point,tetra_P2,tetra_P3,tetra_P4 );
    tetra_V=tetra_volume{ii,1};
    if ((tot_V-tetra_V)<=1e-6)
        p_tetra_num=ii;
        return;
    end
end



end

