function [ p_forward, n_forward, p_backward, n_backward, chk_out ] = intercept_line_tetra_fast( line_point, line_direction, mesh_num, p,t_sort, nei_list_coplane)
%Weizong Xu, August, 2017


%tic
n=0;
p_out=cell(2,1);
n_out=[0,0];
p_forward=[];
p_backward=[];
n_forward=[];
n_backward=[];

tetra_P1=p(:,t_sort(1,mesh_num))';
tetra_P2=p(:,t_sort(2,mesh_num))';
tetra_P3=p(:,t_sort(3,mesh_num))';
tetra_P4=p(:,t_sort(4,mesh_num))';

%toc

%tic
%P1P2P3
[flag, t] = rayTriangleIntersection_fast2(line_point, line_direction, tetra_P1, tetra_P2, tetra_P3);
if ( flag == 1 || flag == 2) %intercept is within the tetrahdron plane
    n=n+1;
    p_out{n}=line_point + t*line_direction;
    n_out(n)=nei_list_coplane.four_neighbour(mesh_num,1);
end

%P1P2P4
[flag, t] = rayTriangleIntersection_fast2(line_point, line_direction, tetra_P1, tetra_P2, tetra_P4);
if ( flag == 1 || flag == 2) %intercept is within the tetrahdron plane
    n=n+1;
    p_out{n}=line_point + t*line_direction;
    n_out(n)=nei_list_coplane.four_neighbour(mesh_num,2);
end

%P1P3P4
[flag, t] = rayTriangleIntersection_fast2(line_point, line_direction, tetra_P1, tetra_P3, tetra_P4);
if ( flag == 1) %intercept is within the tetrahdron plane
    n=n+1;
    p_out{n}=line_point + t*line_direction;
    n_out(n)=nei_list_coplane.four_neighbour(mesh_num,3);
end

%P2P3P4
[flag, t] = rayTriangleIntersection_fast2(line_point, line_direction, tetra_P2, tetra_P3, tetra_P4);
if ( flag == 1) %intercept is within the tetrahdron plane
    n=n+1;
    p_out{n}=line_point + t*line_direction;
    n_out(n)=nei_list_coplane.four_neighbour(mesh_num,4);
end
%toc

%tic
if (n==2) %normally, only intercept two planes
    chk_out=1; %normal return
    p_chk_d=p_out{1}-p_out{2};
%     if (p_chk_d*line_direction'>0)
    if ((p_chk_d(1)*line_direction(1)+p_chk_d(2)*line_direction(2)+p_chk_d(3)*line_direction(3))>0)
        p_forward=p_out{1};
        p_backward=p_out{2};
        n_forward=n_out(1);
        n_backward=n_out(2);
    else
        p_forward=p_out{2};
        p_backward=p_out{1};
        n_forward=n_out(2);
        n_backward=n_out(1);
    end
    
else

if (n==3) %hit at the two plane intersect
    chk_out=2;
    disp('WARNING! Hit at the two plane intersect')
    if isequal(p_out{1},p_out{2}==0)
       p_chk_d=p_out{1}-p_out{2}; 
        if (p_chk_d*line_direction'>0)
            p_forward=p_out{1};
            p_backward=p_out{2};
        else
            p_forward=p_out{2};
            p_backward=p_out{1};
            n_forward=n_out(2);
            n_backward=n_out(1);
        end
    else
        p_chk_d=p_out{1}-p_out{3};
        if (p_chk_d*line_direction'>0)
            p_forward=p_out{1};
            p_backward=p_out{3};
        else
            p_forward=p_out{3};
            p_backward=p_out{1};
        end
    end
    
end

if (n==4) %hit one corner
    chk_out=3;
    disp('WARNING! Hit at the tetrahedron corner')
    if (isequal(p_out{1},p_out{2})==0)
       p_chk_d=p_out{1}-p_out{2}; 
        if (p_chk_d*l_d'>0)
            p_forward=p_out{1};
            p_backward=p_out{2};
        else
            p_forward=p_out{2};
            p_backward=p_out{1};
        end
    end
        
    if (isequal(p_out{1},p_out{3})==0)
       p_chk_d=p_out{1}-p_out{3}; 
        if (p_chk_d*l_d'>0)
            p_forward=p_out{1};
            p_backward=p_out{3};
        else
            p_forward=p_out{3};
            p_backward=p_out{1};
        end
    end   
        
    if (isequal(p_out{1},p_out{4})==0)
       p_chk_d=p_out{1}-p_out{4}; 
        if (p_chk_d*l_d'>0)
            p_forward=p_out{1};
            p_backward=p_out{4};
        else
            p_forward=p_out{4};
            p_backward=p_out{1};
        end
    end
    %only three condition needed to judge

end

if (n<=1)
    chk_out=0; % no intercept point
    p_forward=[];
    n_forward=[];
    p_backward=[];
    n_backward=[];
end



end

%toc

end

