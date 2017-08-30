function [ ] = look_model_tetra( model_mesh_nei, look_p, add_on_coor )
%Weizong Xu, August, 2017

figure;hold on;
tmp_pp=[];pp=[];
for ii=1:length(look_p)
    nei_list=model_mesh_nei.t(1:4,look_p(ii));
    for i=1:4
        tmp_pp(:,i)=model_mesh_nei.p(:,nei_list(i));
    end
    plot3([tmp_pp(1,1) tmp_pp(1,2)],[tmp_pp(2,1) tmp_pp(2,2)],[tmp_pp(3,1) tmp_pp(3,2)],'b')
    plot3([tmp_pp(1,1) tmp_pp(1,3)],[tmp_pp(2,1) tmp_pp(2,3)],[tmp_pp(3,1) tmp_pp(3,3)],'b')
    plot3([tmp_pp(1,1) tmp_pp(1,4)],[tmp_pp(2,1) tmp_pp(2,4)],[tmp_pp(3,1) tmp_pp(3,4)],'b')     
    plot3([tmp_pp(1,2) tmp_pp(1,3)],[tmp_pp(2,2) tmp_pp(2,3)],[tmp_pp(3,2) tmp_pp(3,3)],'b')
    plot3([tmp_pp(1,2) tmp_pp(1,4)],[tmp_pp(2,2) tmp_pp(2,4)],[tmp_pp(3,2) tmp_pp(3,4)],'b')
    plot3([tmp_pp(1,3) tmp_pp(1,4)],[tmp_pp(2,3) tmp_pp(2,4)],[tmp_pp(3,3) tmp_pp(3,4)],'b')
    pp=[pp tmp_pp];
end
%figure;
scatter3(pp(1,:),pp(2,:),pp(3,:),'b')
%hold on;

grid on;
if (isempty(add_on_coor)==0)
    scatter3(add_on_coor(1,:),add_on_coor(2,:),add_on_coor(3,:),'r')
end
hold off;
end

