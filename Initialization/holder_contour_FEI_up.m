function [ FEI_frame_out] = holder_contour_FEI_up()
% input contour of holder frame, constructure a numeric model of it
%Exp calibrated for FEI LB holder.
%Weizong, April, 2015, wxu4@ncsu.edu;

scale_bar = 0.0172766;
p_level1=-0.150/scale_bar;
p_level2=0.100/scale_bar;
p_level3=0.220/scale_bar;
FEI_frame_input = [...
    417,198,p_level1; ...
    317,198,p_level1; ...
    312,198,p_level1; ...
    308,199,p_level1; ...
    301,202,p_level1; ...
    292,206,p_level1; ...
    285,210,p_level1; ...
    277,215,p_level1; ...
    272,220,p_level1; ...
    267,226,p_level1; ...
    263,232,p_level1; ...
    259,239,p_level1; ...
    255,247,p_level1; ...
    254,252,p_level1; ...
    253,393,p_level1; ...
    254,400,p_level1; ...
    256,406,p_level1; ...
    259,412,p_level1; ...
    262,420,p_level1; ...
    267,426,p_level1; ...
    271,431,p_level1; ...
    277,438,p_level1; ...
    284,443,p_level1; ...
    289,447,p_level1; ...
    298,451,p_level1; ...
    306,454,p_level1; ...
    314,457,p_level1; ...
    412,460,p_level1; ...
    412,460,p_level2; ...
    577,460,p_level2; ...
    592,456,p_level2; ...
    598,452,p_level2; ...
    604,449,p_level2; ...
    611,445,p_level2; ...
    611,445,p_level3; ...
    627,435,p_level3; ...
    757,434,p_level3; ...
    757,226,p_level3; ...
    627,224,p_level3; ...
    613,213,p_level3; ...
    613,213,p_level2; ...
    604,207,p_level2; ...
    595,203,p_level2; ...
    588,200,p_level2; ...
    580,198,p_level2; ...
    417,198,p_level2; ...
    417,198,p_level1];

z_offset=0;
FEI_frame_input(:,2)=-FEI_frame_input(:,2); %flip due to image + is y-

%center_coor= [463,-324,0]; %old version
center_coor= [462,-328,0]; %newly version
FEI_frame_coor = zeros(length(FEI_frame_input),3);
FEI_frame_coor(:,1) = (FEI_frame_input(:,1)-center_coor(1))*scale_bar;
FEI_frame_coor(:,2) = (FEI_frame_input(:,2)-center_coor(2))*scale_bar;
FEI_frame_coor(:,3) = (FEI_frame_input(:,3)-center_coor(3))*scale_bar+z_offset;
%figure;plot3(FEI_frame_coor(:,1),FEI_frame_coor(:,2),FEI_frame_coor(:,3))
%figure;scatter3(FEI_frame_coor(:,1),FEI_frame_coor(:,2),FEI_frame_coor(:,3))
%figure;scatter3(FEI_frame_input(:,1),FEI_frame_input(:,2),FEI_frame_input(:,3))

d_step_xy_max=0.1;
d_step_z_max =0.5;
chk=2;
FEI_frame_coor2=FEI_frame_coor;
while (chk==2)

    clear FEI_frame_out;
    num=0;
    chk=1;
for i=1:length(FEI_frame_coor2)-1
    num=num+1;
    FEI_frame_out(num,1)=FEI_frame_coor2(i,1); %input first one
    FEI_frame_out(num,2)=FEI_frame_coor2(i,2);
    FEI_frame_out(num,3)=FEI_frame_coor2(i,3);

    %intercept
    dist_chk=sqrt((FEI_frame_coor2(i,1)-FEI_frame_coor2(i+1,1))^2 + ...
        (FEI_frame_coor2(i,2)-FEI_frame_coor2(i+1,2))^2+ ...
        (FEI_frame_coor2(i,3)-FEI_frame_coor2(i+1,3))^2);
    z_chk=abs(FEI_frame_coor2(i,3)-FEI_frame_coor2(i+1,3));
    if (((dist_chk > d_step_xy_max) && (z_chk==0)) || ((dist_chk > d_step_z_max) && (z_chk>0)))
        num=num+1;
        FEI_frame_out(num,1)=(FEI_frame_coor2(i,1)+FEI_frame_coor2(i+1,1))/2;
        FEI_frame_out(num,2)=(FEI_frame_coor2(i,2)+FEI_frame_coor2(i+1,2))/2;
        FEI_frame_out(num,3)=(FEI_frame_coor2(i,3)+FEI_frame_coor2(i+1,3))/2;
        chk=2;
    end

end
        FEI_frame_out(num+1,1)=FEI_frame_coor2(i+1,1); %input last one
        FEI_frame_out(num+1,2)=FEI_frame_coor2(i+1,2);
        FEI_frame_out(num+1,3)=FEI_frame_coor2(i+1,3);
        FEI_frame_coor2=FEI_frame_out;
end

%get plane normal
temp=length(FEI_frame_out);
FEI_frame_out(:,4)=zeros(temp,1);
FEI_frame_out(:,5)=zeros(temp,1);
FEI_frame_out(:,6)=zeros(temp,1);

for i=1:temp-1
    V_x=FEI_frame_out(i,1)-FEI_frame_out(i+1,1);
    V_y=FEI_frame_out(i,2)-FEI_frame_out(i+1,2);
    V_z=FEI_frame_out(i,3)-FEI_frame_out(i+1,3);
    if (V_z==0)
        V_cross=cross([V_x,V_y,0],[0,0,1]);
        V_cross=V_cross/norm(V_cross);
        FEI_frame_out(i,4)=V_cross(1);
        FEI_frame_out(i,5)=V_cross(2);
        FEI_frame_out(i,6)=V_cross(3);
    else
        FEI_frame_out(i,4)=FEI_frame_out(i-1,4);
        FEI_frame_out(i,5)=FEI_frame_out(i-1,5);
        FEI_frame_out(i,6)=FEI_frame_out(i-1,6);
    end
    %last point must be considered
    FEI_frame_out(temp,4)=FEI_frame_out(temp-1,4);
    FEI_frame_out(temp,5)=FEI_frame_out(temp-1,5);
    FEI_frame_out(temp,6)=FEI_frame_out(temp-1,6);
end

%special points should be considered
 temp=length(FEI_frame_out);
 FEI_frame_out(:,7)=ones(temp,1);

for i=1:temp-1
    if (abs(FEI_frame_out(i,3)-FEI_frame_out(i+1,3))>0.0001 && FEI_frame_out(i,1)<0)
        if (FEI_frame_out(i,3)>FEI_frame_out(i+1,3)) %only front surface point is considered
            FEI_frame_out(i,7)=2; % =2 analyse using vertical plane
        else
            FEI_frame_out(i+1,7)=2;
        end
    end
    %=1 others use horizonal plane
end

for i=1:temp-1
    if (abs(FEI_frame_out(i,3)-FEI_frame_out(i+1,3))>0.0001 && FEI_frame_out(i,1)>0 && FEI_frame_out(i+1,1)>0)
        if (FEI_frame_out(i,3)<FEI_frame_out(i+1,3)) %only front surface point is considered
            FEI_frame_out(i,7)=3; % =2 analyse using vertical plane
        else
            FEI_frame_out(i+1,7)=3;
        end
    end
    %=1 others use horizonal plane
end
%*****************************




%for display only
d_step_xy_max=0.1;
d_step_z_max =0.5;
chk=2;
FEI_frame_coor2=FEI_frame_coor;
while (chk==2)

    clear FEI_frame_disp;
    num=0;
    chk=1;
for i=1:length(FEI_frame_coor2)-1
    num=num+1;
    FEI_frame_disp(num,1)=FEI_frame_coor2(i,1); %input first one
    FEI_frame_disp(num,2)=FEI_frame_coor2(i,2);
    FEI_frame_disp(num,3)=FEI_frame_coor2(i,3);

    %intercept
    dist_chk=sqrt((FEI_frame_coor2(i,1)-FEI_frame_coor2(i+1,1))^2 + ...
        (FEI_frame_coor2(i,2)-FEI_frame_coor2(i+1,2))^2+ ...
        (FEI_frame_coor2(i,3)-FEI_frame_coor2(i+1,3))^2);
    z_chk=abs(FEI_frame_coor2(i,3)-FEI_frame_coor2(i+1,3));
    if (((dist_chk > d_step_xy_max) && (z_chk==0)) || ((dist_chk > d_step_z_max) && (z_chk>0)))
        num=num+1;
        FEI_frame_disp(num,1)=(FEI_frame_coor2(i,1)+FEI_frame_coor2(i+1,1))/2;
        FEI_frame_disp(num,2)=(FEI_frame_coor2(i,2)+FEI_frame_coor2(i+1,2))/2;
        FEI_frame_disp(num,3)=(FEI_frame_coor2(i,3)+FEI_frame_coor2(i+1,3))/2;
        chk=2;
    end

end
        FEI_frame_disp(num+1,1)=FEI_frame_coor2(i+1,1); %input last one
        FEI_frame_disp(num+1,2)=FEI_frame_coor2(i+1,2);
        FEI_frame_disp(num+1,3)=FEI_frame_coor2(i+1,3);
        FEI_frame_coor2=FEI_frame_disp;
end

%topo to 3D
temp_length = length (FEI_frame_disp);
z_step=0.02;
num=0;
for i=1:temp_length
    num=num+1;
    FEI_frame_disp_3D(num,1)=FEI_frame_disp(i,1);
    FEI_frame_disp_3D(num,2)=FEI_frame_disp(i,2);
    FEI_frame_disp_3D(num,3)=FEI_frame_disp(i,3);
    t_chk=FEI_frame_disp(i,3);
    while (t_chk-z_step>-0.3)
        num=num+1;
        FEI_frame_disp_3D(num,1)=FEI_frame_disp(i,1);
        FEI_frame_disp_3D(num,2)=FEI_frame_disp(i,2);
        FEI_frame_disp_3D(num,3)=FEI_frame_disp_3D(num-1,3)-z_step;
        t_chk=FEI_frame_disp_3D(num,3);
    end
end

% Plane1_P1 = [417,198,p_level2];
% Plane1_P2 = [417,198,-p_level2*10];
% Plane1_P3 = [417,140-50,-p_level2*10];
% Plane1_P4 = [417,140-50,p_level2];
Plane1_P1 = [417,198,p_level2];
Plane1_P2 = [417,140-50,p_level2];
Plane1_P3 = [417,140-50,p_level1-1.0/scale_bar];
Plane1_P4 = [417,198,p_level1-1.0/scale_bar];
Plane1 = [Plane1_P1;Plane1_P2;Plane1_P3;Plane1_P4;Plane1_P1]; %last one is added for display
Plane1(:,2)=-Plane1(:,2); %flip due to image + is y-

Plane1_coor = zeros(5,3);
Plane1_coor(:,1) = (Plane1(:,1)-center_coor(1))*scale_bar;
Plane1_coor(:,2) = (Plane1(:,2)-center_coor(2))*scale_bar;
Plane1_coor(:,3) = (Plane1(:,3)-center_coor(3))*scale_bar+z_offset;

Plane2_P1 = [613,213,p_level3];
Plane2_P2 = [613,213,p_level2-1.0/scale_bar];
Plane2_P3 = [613,165,p_level2-1.0/scale_bar];
Plane2_P4 = [613,165,p_level3];
Plane2 = [Plane2_P1;Plane2_P2;Plane2_P3;Plane2_P4;Plane2_P1];
Plane2(:,2)=-Plane2(:,2); %flip due to image + is y-

Plane2_coor = zeros(5,3);
Plane2_coor(:,1) = (Plane2(:,1)-center_coor(1))*scale_bar;
Plane2_coor(:,2) = (Plane2(:,2)-center_coor(2))*scale_bar;
Plane2_coor(:,3) = (Plane2(:,3)-center_coor(3))*scale_bar+z_offset;

Plane3_P1 = [611,445,p_level3];
Plane3_P2 = [611,489,p_level3];
Plane3_P3 = [611,489,p_level2-1.0/scale_bar];
Plane3_P4 = [611,445,p_level2-1.0/scale_bar];
Plane3 = [Plane3_P1;Plane3_P2;Plane3_P3;Plane3_P4;Plane3_P1];
Plane3(:,2)=-Plane3(:,2); %flip due to image + is y-

Plane3_coor = zeros(5,3);
Plane3_coor(:,1) = (Plane3(:,1)-center_coor(1))*scale_bar;
Plane3_coor(:,2) = (Plane3(:,2)-center_coor(2))*scale_bar;
Plane3_coor(:,3) = (Plane3(:,3)-center_coor(3))*scale_bar+z_offset;

Plane4_P1 = [412,460,p_level2];
Plane4_P2 = [412,460,p_level1-1.0/scale_bar];
Plane4_P3 = [412,512+50,p_level1-1.0/scale_bar];
Plane4_P4 = [412,512+50,p_level2];
Plane4 = [Plane4_P1;Plane4_P2;Plane4_P3;Plane4_P4;Plane4_P1];
Plane4(:,2)=-Plane4(:,2); %flip due to image + is y-

Plane4_coor = zeros(5,3);
Plane4_coor(:,1) = (Plane4(:,1)-center_coor(1))*scale_bar;
Plane4_coor(:,2) = (Plane4(:,2)-center_coor(2))*scale_bar;
Plane4_coor(:,3) = (Plane4(:,3)-center_coor(3))*scale_bar+z_offset;


Plane_coor(:,:,1)=Plane1_coor;
Plane_coor(:,:,2)=Plane2_coor;
Plane_coor(:,:,3)=Plane3_coor;
Plane_coor(:,:,4)=Plane4_coor;

temp=size(Plane_coor);
num=1;
for i=1:temp(3) %Pack these plane normal into one output matrix
    for j=1:4
        for k=1:3
            FEI_frame_out(num,8)=Plane_coor(j,k,i);
            num=num+1;
        end
    end
end

% figure;hold on;
% plot3(Plane1_coor(:,1),Plane1_coor(:,2),Plane1_coor(:,3),'b')
% plot3(Plane2_coor(:,1),Plane2_coor(:,2),Plane2_coor(:,3),'b')
% plot3(Plane3_coor(:,1),Plane3_coor(:,2),Plane3_coor(:,3),'b')
% plot3(Plane4_coor(:,1),Plane4_coor(:,2),Plane4_coor(:,3),'b')
% scatter3(0,0,0,200,'filled','d','r');
% scatter3(FEI_frame_out(:,1),FEI_frame_out(:,2),FEI_frame_out(:,3),'b')
% grid on
% hold off;




%figure;scatter3(FEI_frame_out(:,1),FEI_frame_out(:,2),FEI_frame_out(:,3))
%figure;scatter3(FEI_frame_disp(:,1),FEI_frame_disp(:,2),FEI_frame_disp(:,3))
%figure;scatter3(FEI_frame_disp_3D(:,1),FEI_frame_disp_3D(:,2),FEI_frame_disp_3D(:,3))
end

