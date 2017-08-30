function [ p_forward, n_forward, p_backward, n_backward, chk_out ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_mesh_nei )
%Weizong Xu, August, 2017


%tic
p_out=[];n_out=[];p_forward=[];p_backward=[];n_forward=[];n_backward=[];
tetra_P1=model_mesh_nei.p(:,model_mesh_nei.t_sort(1,mesh_num))';
tetra_P2=model_mesh_nei.p(:,model_mesh_nei.t_sort(2,mesh_num))';
tetra_P3=model_mesh_nei.p(:,model_mesh_nei.t_sort(3,mesh_num))';
tetra_P4=model_mesh_nei.p(:,model_mesh_nei.t_sort(4,mesh_num))';

nei_P123=model_mesh_nei.nei_list_coplane.four_neighbour(mesh_num,1);
nei_P124=model_mesh_nei.nei_list_coplane.four_neighbour(mesh_num,2);
nei_P134=model_mesh_nei.nei_list_coplane.four_neighbour(mesh_num,3);
nei_P234=model_mesh_nei.nei_list_coplane.four_neighbour(mesh_num,4);
%toc
%tic
l_d=line_direction/norm(line_direction); %line direction %u
p_n123=cross(tetra_P1-tetra_P2, tetra_P1-tetra_P3);
p_n124=cross(tetra_P1-tetra_P2, tetra_P1-tetra_P4);
p_n134=cross(tetra_P1-tetra_P3, tetra_P1-tetra_P4);
p_n234=cross(tetra_P2-tetra_P3, tetra_P2-tetra_P4);
%toc
%tic
l_p123=line_point-tetra_P1; %w
if (l_p123==0) l_p123=line_point-tetra_P2; end %if the start point happen to be the same point in the tetrahedron corner
%P123 = line_point-dot(p_n123,l_p123)/dot(p_n123,l_d).*l_d;
P123 = line_point-(p_n123*l_p123')/(p_n123*l_d').*l_d;
%if if this point is within the tetrahedron surface
v0=tetra_P3-tetra_P1;
v1=tetra_P2-tetra_P1;
v2=P123-tetra_P1;
dot00 = v0*v0';%dot(v0, v0);
dot01 = v0*v1';%dot(v0, v1);
dot02 = v0*v2';%dot(v0, v2);
dot11 = v1*v1';%dot(v1, v1);
dot12 = v1*v2';%dot(v1, v2);
invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
u = (dot11 * dot02 - dot01 * dot12) * invDenom;
v = (dot00 * dot12 - dot01 * dot02) * invDenom;
if (u >= 0 && v >= 0 && u + v <= 1) %intercept is within the tetrahdron plane
    p_out=[p_out;P123];
    n_out=[n_out;nei_P123];
end
%toc
%tic
l_p124=line_point-tetra_P1; %w
if (l_p124==0) l_p124=line_point-tetra_P2; end %if the start point happen to be the same point in the tetrahedron corner
%P124 = line_point-dot(p_n124,l_p124)/dot(p_n124,l_d).*l_d;
P124 = line_point-(p_n124*l_p124')/(p_n124*l_d').*l_d;
%if if this point is within the tetrahedron surface
v0=tetra_P4-tetra_P1;
v1=tetra_P2-tetra_P1;
v2=P124-tetra_P1;
dot00 = v0*v0';%dot(v0, v0);
dot01 = v0*v1';%dot(v0, v1);
dot02 = v0*v2';%dot(v0, v2);
dot11 = v1*v1';%dot(v1, v1);
dot12 = v1*v2';%dot(v1, v2);
invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
u = (dot11 * dot02 - dot01 * dot12) * invDenom;
v = (dot00 * dot12 - dot01 * dot02) * invDenom;
if (u >= 0 && v >= 0 && u + v <= 1) %intercept is within the tetrahdron plane
    p_out=[p_out;P124];
    n_out=[n_out;nei_P124];
end
%toc
%tic
l_p134=line_point-tetra_P1; %w
if (l_p134==0) l_p134=line_point-tetra_P3; end %if the start point happen to be the same point in the tetrahedron corner
%P134 = line_point-dot(p_n134,l_p134)/dot(p_n134,l_d).*l_d;
P134 = line_point-(p_n134*l_p134')/(p_n134*l_d').*l_d;
%if if this point is within the tetrahedron surface
v0=tetra_P4-tetra_P1;
v1=tetra_P3-tetra_P1;
v2=P134-tetra_P1;
dot00 = v0*v0';%dot(v0, v0);
dot01 = v0*v1';%dot(v0, v1);
dot02 = v0*v2';%dot(v0, v2);
dot11 = v1*v1';%dot(v1, v1);
dot12 = v1*v2';%dot(v1, v2);
invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
u = (dot11 * dot02 - dot01 * dot12) * invDenom;
v = (dot00 * dot12 - dot01 * dot02) * invDenom;
if (u >= 0 && v >= 0 && u + v <= 1) %intercept is within the tetrahdron plane
    p_out=[p_out;P134];
    n_out=[n_out;nei_P134];
end
%toc
%%tic
l_p234=line_point-tetra_P2; %w
if (l_p234==0) l_p234=line_point-tetra_P3; end %if the start point happen to be the same point in the tetrahedron corner
%P234 = line_point-dot(p_n234,l_p234)/dot(p_n234,l_d).*l_d;
P234 = line_point-(p_n234*l_p234')/(p_n234*l_d').*l_d;
%if if this point is within the tetrahedron surface
v0=tetra_P4-tetra_P2;
v1=tetra_P3-tetra_P2;
v2=P234-tetra_P2;
dot00 = v0*v0';%dot(v0, v0);
dot01 = v0*v1';%dot(v0, v1);
dot02 = v0*v2';%dot(v0, v2);
dot11 = v1*v1';%dot(v1, v1);
dot12 = v1*v2';%dot(v1, v2);
invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
u = (dot11 * dot02 - dot01 * dot12) * invDenom;
v = (dot00 * dot12 - dot01 * dot02) * invDenom;
if (u >= 0 && v >= 0 && u + v <= 1) %intercept is within the tetrahdron plane
    p_out=[p_out;P234];
    n_out=[n_out;nei_P234];
end
%toc
%tic
if (size(p_out,1)==2) %normally, only intercept two planes
    chk_out=1; %normal return
    p_chk_d=p_out(1,:)-p_out(2,:);
    if (dot(p_chk_d,l_d)>0)
        p_forward=p_out(1,:);
        p_backward=p_out(2,:);
        n_forward=n_out(1,:);
        n_backward=n_out(2,:);
    else
        p_forward=p_out(2,:);
        p_backward=p_out(1,:);
        n_forward=n_out(2,:);
        n_backward=n_out(1,:);
    end
end

if (size(p_out,1)==3) %hit at the two plane intersect
    chk_out=2;
    disp('WARNING! Hit at the two plane intersect')
    if isequal(p_out(1,:),p_out(2,:)==0)
       p_chk_d=p_out(1,:)-p_out(2,:); 
        if (dot(p_chk_d,l_d)>0)
            p_forward=p_out(1,:);
            p_backward=p_out(2,:);
        else
            p_forward=p_out(2,:);
            p_backward=p_out(1,:);
        end
    else
        p_chk_d=p_out(1,:)-p_out(3,:);
        if (dot(p_chk_d,l_d)>0)
            p_forward=p_out(1,:);
            p_backward=p_out(3,:);
        else
            p_forward=p_out(3,:);
            p_backward=p_out(1,:);
        end
    end
end


if (size(p_out,1)==4) %hit one corner
    chk_out=3;
    disp('WARNING! Hit at the tetrahedron corner')
    if (isequal(p_out(1,:),p_out(2,:))==0)
       p_chk_d=p_out(1,:)-p_out(2,:); 
        if (dot(p_chk_d,l_d)>0)
            p_forward=p_out(1,:);
            p_backward=p_out(2,:);
        else
            p_forward=p_out(2,:);
            p_backward=p_out(1,:);
        end
    end
        
    if (isequal(p_out(1,:),p_out(3,:))==0)
       p_chk_d=p_out(1,:)-p_out(3,:); 
        if (dot(p_chk_d,l_d)>0)
            p_forward=p_out(1,:);
            p_backward=p_out(3,:);
        else
            p_forward=p_out(3,:);
            p_backward=p_out(1,:);
        end
    end   
        
    if (isequal(p_out(1,:),p_out(4,:))==0)
       p_chk_d=p_out(1,:)-p_out(4,:); 
        if (dot(p_chk_d,l_d)>0)
            p_forward=p_out(1,:);
            p_backward=p_out(4,:);
        else
            p_forward=p_out(4,:);
            p_backward=p_out(1,:);
        end
    end
    %only three condition needed to judge

end

if (size(p_out,1)<=1)
    chk_out=0; % no intercept point
    p_forward=[];
    n_forward=[];
    p_backward=[];
    n_backward=[];
end
%Many thanks
% % %source http://www.mathworks.com/matlabcentral/fileexchange/17751-straight-line-and-plane-intersection
% % % n: normal vector of the Plane 
% % % V0: any point that belongs to the Plane 
% % % P0: end point 1 of the segment P0P1
% % % P1:  end point 2 of the segment P0P1
% % % I=[0 0 0];
% % % u = P1-P0;
% % % w = P0 - V0;
% % % D = dot(n,u);
% % % N = -dot(n,w);
% % % check=0;
% % % if abs(D) < 10^-7        % The segment is parallel to plane
% % %         if N == 0           % The segment lies in plane
% % %             check=2;
% % %             return
% % %         else
% % %             check=0;       %no intersection
% % %             return
% % %         end
% % % end
% % % 
% % % %compute the intersection parameter
% % % sI = N / D;
% % % I = P0+ sI.*u;
% % % 
% % % if (sI < 0 || sI > 1)
% % %     check= 3;          %The intersection point  lies outside the segment, so there is no intersection
% % % else
% % %     check=1;
% % % end
% % %P123 = line_point-dot(p_n123,l_p)/dot(p_n123,l_d).*l_d;
% % 
% % %http://www.blackpawn.com/texts/pointinpoly/
% % % %%Compute vectors        
% % % v0 = C - A
% % % v1 = B - A
% % % v2 = P - A
% % % 
% % % %Compute dot products
% % % dot00 = dot(v0, v0)
% % % dot01 = dot(v0, v1)
% % % dot02 = dot(v0, v2)
% % % dot11 = dot(v1, v1)
% % % dot12 = dot(v1, v2)
% % % 
% % % %Compute barycentric coordinates
% % % invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
% % % u = (dot11 * dot02 - dot01 * dot12) * invDenom
% % % v = (dot00 * dot12 - dot01 * dot02) * invDenom
% % % 
% % % %Check if point is in triangle
% % % return (u >= 0) && (v >= 0) && (u + v < 1)


%toc
%look_model_tetra( model_mesh_nei, [mesh_num n_forward], [line_point;p_out]' );
end

