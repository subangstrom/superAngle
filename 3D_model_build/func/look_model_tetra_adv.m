function [ ] = look_model_tetra_adv( model_mesh_nei, look_p, add_on_coor, look_p2, add_on_coor2, x_range, y_range,ether_start_p )
%Weizong Xu, August, 2017
%   Detailed explanation goes here
%color1=[82 152 197]/255;
%color12=[0 152 158]/255;
%color12=[15 198 92]/255;
% color1=[140 198 98]/255;
color1=[3 110 184]/255;
color3=[233 56 5]/255;
% color2=[255 128 0]/255;
color2=[216 34 13]/255;
color12=[246 172 26]/255;
%color2=[0.5,0.5,0.5];%gray
%**************************************Draw background*******************************************
look_p_all=1:length(model_mesh_nei.t_sort);
look_p_all_reduced=[];
for i=1:length(look_p_all)
    %look_p_all_all(i)->t_sort
    t_i=model_mesh_nei.t_sort(:,look_p_all(i));
    chk_add=1;
    for j=1:4
        p_j=model_mesh_nei.p(:,t_i(j));
        if (p_j(1)<x_range(1) || p_j(1)>x_range(2) || p_j(2)<y_range(1) || p_j(2)>y_range(2))
            chk_add=0;
        end
    end
    if (chk_add==1)
        look_p_all_reduced=[look_p_all_reduced look_p_all(i)];
    end
end
figure;hold on;
tmp_pp=[];
for ii=1:length(look_p_all_reduced)
    nei_list=model_mesh_nei.t(1:4,look_p_all_reduced(ii));
    for i=1:4
        tmp_pp(:,i)=model_mesh_nei.p(:,nei_list(i));
    end
    factor=1.5;
    plot3([tmp_pp(1,1) tmp_pp(1,2)],[tmp_pp(2,1) tmp_pp(2,2)],[tmp_pp(3,1) tmp_pp(3,2)],'color', [0.5,0.5,0.5]*factor)
    plot3([tmp_pp(1,1) tmp_pp(1,3)],[tmp_pp(2,1) tmp_pp(2,3)],[tmp_pp(3,1) tmp_pp(3,3)],'color', [0.5,0.5,0.5]*factor)
    plot3([tmp_pp(1,1) tmp_pp(1,4)],[tmp_pp(2,1) tmp_pp(2,4)],[tmp_pp(3,1) tmp_pp(3,4)],'color', [0.5,0.5,0.5]*factor)     
    plot3([tmp_pp(1,2) tmp_pp(1,3)],[tmp_pp(2,2) tmp_pp(2,3)],[tmp_pp(3,2) tmp_pp(3,3)],'color', [0.5,0.5,0.5]*factor)
    plot3([tmp_pp(1,2) tmp_pp(1,4)],[tmp_pp(2,2) tmp_pp(2,4)],[tmp_pp(3,2) tmp_pp(3,4)],'color', [0.5,0.5,0.5]*factor)
    plot3([tmp_pp(1,3) tmp_pp(1,4)],[tmp_pp(2,3) tmp_pp(2,4)],[tmp_pp(3,3) tmp_pp(3,4)],'color', [0.5,0.5,0.5]*factor)
end

%************************************Draw slice model ************************************************
tmp_pp=[];pp=[];
for ii=1:length(look_p)
    nei_list=model_mesh_nei.t(1:4,look_p(ii));
    for i=1:4
        tmp_pp(:,i)=model_mesh_nei.p(:,nei_list(i));
    end
    plot3([tmp_pp(1,1) tmp_pp(1,2)],[tmp_pp(2,1) tmp_pp(2,2)],[tmp_pp(3,1) tmp_pp(3,2)],'color',color2,'LineWidth',0.5)
    plot3([tmp_pp(1,1) tmp_pp(1,3)],[tmp_pp(2,1) tmp_pp(2,3)],[tmp_pp(3,1) tmp_pp(3,3)],'color',color2,'LineWidth',0.5)
    plot3([tmp_pp(1,1) tmp_pp(1,4)],[tmp_pp(2,1) tmp_pp(2,4)],[tmp_pp(3,1) tmp_pp(3,4)],'color',color2,'LineWidth',0.5)     
    plot3([tmp_pp(1,2) tmp_pp(1,3)],[tmp_pp(2,2) tmp_pp(2,3)],[tmp_pp(3,2) tmp_pp(3,3)],'color',color2,'LineWidth',0.5)
    plot3([tmp_pp(1,2) tmp_pp(1,4)],[tmp_pp(2,2) tmp_pp(2,4)],[tmp_pp(3,2) tmp_pp(3,4)],'color',color2,'LineWidth',0.5)
    plot3([tmp_pp(1,3) tmp_pp(1,4)],[tmp_pp(2,3) tmp_pp(2,4)],[tmp_pp(3,3) tmp_pp(3,4)],'color',color2,'LineWidth',0.5)
    pp=[pp tmp_pp];
end
%scatter3(pp(1,:),pp(2,:),pp(3,:),'g')
% scatter3(pp(1,:),pp(2,:),pp(3,:),'MarkerEdgeColor',color2,...
%               'MarkerFaceColor',color2,...
%               'LineWidth',1.5)
% grid on;
if (~isempty(add_on_coor))
    %scatter3(add_on_coor(1,:),add_on_coor(2,:),add_on_coor(3,:),'.r')
    scatter3(add_on_coor(1,:),add_on_coor(2,:),add_on_coor(3,:),10,'filled','MarkerFaceColor',color3)
end

%************************************Draw beam ************************************************
line1=ether_start_p(1);
line2=ether_start_p(2);
line30=ether_start_p(3);
line31=ether_start_p(3)-10+6;
line32=ether_start_p(3)+20-6;
plot3([line1,line1],[line2,line2],[line31,line30],'LineWidth',2,'color',color3)
plot3([line1,line1],[line2,line2],[line32,line30+10],'LineWidth',2,'color',color3)

%************************************Draw exit X-ray ************************************************
if ~isempty(look_p2)

tmp_pp=[];pp=[];
for ii=1:length(look_p2)
    nei_list=model_mesh_nei.t(1:4,look_p2(ii));
    for i=1:4
        tmp_pp(:,i)=model_mesh_nei.p(:,nei_list(i));
    end
    plot3([tmp_pp(1,1) tmp_pp(1,2)],[tmp_pp(2,1) tmp_pp(2,2)],[tmp_pp(3,1) tmp_pp(3,2)],'color',color1,'LineWidth',0.75)
    plot3([tmp_pp(1,1) tmp_pp(1,3)],[tmp_pp(2,1) tmp_pp(2,3)],[tmp_pp(3,1) tmp_pp(3,3)],'color',color1,'LineWidth',0.75)
    plot3([tmp_pp(1,1) tmp_pp(1,4)],[tmp_pp(2,1) tmp_pp(2,4)],[tmp_pp(3,1) tmp_pp(3,4)],'color',color1,'LineWidth',0.75)     
    plot3([tmp_pp(1,2) tmp_pp(1,3)],[tmp_pp(2,2) tmp_pp(2,3)],[tmp_pp(3,2) tmp_pp(3,3)],'color',color1,'LineWidth',0.75)
    plot3([tmp_pp(1,2) tmp_pp(1,4)],[tmp_pp(2,2) tmp_pp(2,4)],[tmp_pp(3,2) tmp_pp(3,4)],'color',color1,'LineWidth',0.75)
    plot3([tmp_pp(1,3) tmp_pp(1,4)],[tmp_pp(2,3) tmp_pp(2,4)],[tmp_pp(3,3) tmp_pp(3,4)],'color',color1,'LineWidth',0.75)
    pp=[pp tmp_pp];
end
scatter3(pp(1,:),pp(2,:),pp(3,:),'b')
scatter3(pp(1,:),pp(2,:),pp(3,:),'MarkerEdgeColor',color1,...
              'MarkerFaceColor',color1,...
              'LineWidth',1.5)
%grid on;


if (isempty(add_on_coor2)==0)
    %scatter3(add_on_coor2(1,:),add_on_coor2(2,:),add_on_coor2(3,:),'.r')
    plot3(add_on_coor2(1,:),add_on_coor2(2,:),add_on_coor2(3,:),'color',color12,'LineWidth',2)
end

end

zlim([-10 20])
az = -2;
el = 3;
view(az, el);



hold off;
%print(['model_illustration.eps'],'-depsc','-tiff');
%print(['model_illustration'],'-dpng','-r1200');
end

