function [ p_forward, n_forward ] = intercept_line_tetra_fast_forward( line_point, line_direction, tetra_P1,tetra_P2,tetra_P3,tetra_P4, four_neighbour)
%Weizong Xu, August, 2017

%speed up, note: very limited functionality, may have problem when line hit
%corner or edge of tetrahderal

%tic
 % initializing
% tetra_P1=p(:,t_sort(1,mesh_num))';
% tetra_P2=p(:,t_sort(2,mesh_num))';
% tetra_P3=p(:,t_sort(3,mesh_num))';
% tetra_P4=p(:,t_sort(4,mesh_num))';

%toc
% tetra_P1=pp{t_sort(1,mesh_num)};
% tetra_P2=pp{t_sort(2,mesh_num)};
% tetra_P3=pp{t_sort(3,mesh_num)};
% tetra_P4=pp{t_sort(4,mesh_num)};


%P1P2P3
[flag, t] = rayTriangleIntersection_fast2(line_point, line_direction, tetra_P1, tetra_P2, tetra_P3);
if (flag == 1 && t>1e-10) %intercept is within the tetrahdron plane, t>0 foward direction
    p_forward=line_point + t*line_direction;
    %n_forward=four_neighbour(mesh_num,1);
    n_forward=four_neighbour(1);
    return; 
end

%P1P2P4
[flag, t] = rayTriangleIntersection_fast2(line_point, line_direction, tetra_P1, tetra_P2, tetra_P4);
if (flag == 1 && t>1e-10) %intercept is within the tetrahdron plane
    p_forward=line_point + t*line_direction;
    n_forward=four_neighbour(2);
    return;
end

%P1P3P4
[flag, t] = rayTriangleIntersection_fast2(line_point, line_direction, tetra_P1, tetra_P3, tetra_P4);
if (flag == 1 && t>1e-10) %intercept is within the tetrahdron plane
    p_forward=line_point + t*line_direction;
    n_forward=four_neighbour(3);
    return;   
end

%P2P3P4
[flag, t] = rayTriangleIntersection_fast2(line_point, line_direction, tetra_P2, tetra_P3, tetra_P4);
if (flag == 1 && t>1e-10) %intercept is within the tetrahdron plane
    p_forward=line_point + t*line_direction;
    n_forward=four_neighbour(4);
    return;
end
%toc

%return with no intercept point or intercept at corner or line edge
p_forward=[];
n_forward=nan;
end

