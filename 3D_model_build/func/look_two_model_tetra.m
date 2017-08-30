function [ ] = look_two_model_tetra( model_mesh_neiA, look_pA, model_mesh_neiB, look_pB, add_on_coor )
%Weizong Xu, August, 2017

figure;hold on;
tmp_pp=[];pp=[];
for ii=1:length(look_pA)
    nei_list=model_mesh_neiA.t(1:4,look_pA(ii));
    for i=1:4
        tmp_pp(:,i)=model_mesh_neiA.p(:,nei_list(i));
    end
    plot3([tmp_pp(1,1) tmp_pp(1,2)],[tmp_pp(2,1) tmp_pp(2,2)],[tmp_pp(3,1) tmp_pp(3,2)],'b')
    plot3([tmp_pp(1,1) tmp_pp(1,3)],[tmp_pp(2,1) tmp_pp(2,3)],[tmp_pp(3,1) tmp_pp(3,3)],'b')
    plot3([tmp_pp(1,1) tmp_pp(1,4)],[tmp_pp(2,1) tmp_pp(2,4)],[tmp_pp(3,1) tmp_pp(3,4)],'b')     
    plot3([tmp_pp(1,2) tmp_pp(1,3)],[tmp_pp(2,2) tmp_pp(2,3)],[tmp_pp(3,2) tmp_pp(3,3)],'b')
    plot3([tmp_pp(1,2) tmp_pp(1,4)],[tmp_pp(2,2) tmp_pp(2,4)],[tmp_pp(3,2) tmp_pp(3,4)],'b')
    plot3([tmp_pp(1,3) tmp_pp(1,4)],[tmp_pp(2,3) tmp_pp(2,4)],[tmp_pp(3,3) tmp_pp(3,4)],'b')
    pp=[pp tmp_pp];
end
scatter3(pp(1,:),pp(2,:),pp(3,:),'b')

if (isempty(look_pB)==0)
tmp_pp=[];pp=[];
for ii=1:length(look_pB)
    nei_list=model_mesh_neiB.t(1:4,look_pB(ii));
    for i=1:4
        tmp_pp(:,i)=model_mesh_neiB.p(:,nei_list(i));
    end
    plot3([tmp_pp(1,1) tmp_pp(1,2)],[tmp_pp(2,1) tmp_pp(2,2)],[tmp_pp(3,1) tmp_pp(3,2)],'k')
    plot3([tmp_pp(1,1) tmp_pp(1,3)],[tmp_pp(2,1) tmp_pp(2,3)],[tmp_pp(3,1) tmp_pp(3,3)],'k')
    plot3([tmp_pp(1,1) tmp_pp(1,4)],[tmp_pp(2,1) tmp_pp(2,4)],[tmp_pp(3,1) tmp_pp(3,4)],'k')     
    plot3([tmp_pp(1,2) tmp_pp(1,3)],[tmp_pp(2,2) tmp_pp(2,3)],[tmp_pp(3,2) tmp_pp(3,3)],'k')
    plot3([tmp_pp(1,2) tmp_pp(1,4)],[tmp_pp(2,2) tmp_pp(2,4)],[tmp_pp(3,2) tmp_pp(3,4)],'k')
    plot3([tmp_pp(1,3) tmp_pp(1,4)],[tmp_pp(2,3) tmp_pp(2,4)],[tmp_pp(3,3) tmp_pp(3,4)],'k')
    pp=[pp tmp_pp];
end
scatter3(pp(1,:),pp(2,:),pp(3,:),'k')

end

grid on;
if (isempty(add_on_coor)==0)
    scatter3(add_on_coor(1,:),add_on_coor(2,:),add_on_coor(3,:),'r')
end
hold off;
end

