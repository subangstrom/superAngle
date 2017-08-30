function [ Plane1_coor, Plane2_coor ] = holder_carrier_FEI_up()
%Input geometry of FEI LB holder carrier geometry
%input contour of holder frame, constructure a numeric model of it
%Weizong Xu, Apr 2015, wxu4@ncsu.edu;
scale_bar = 0.0172766;
p_level1=-0.22/scale_bar;
p_level2=0;
center_coor= [462,-328,0]; %newly version

Plane1_P1 = [412,244,p_level2];
Plane1_P2 = [412,215,p_level2];
Plane1_P3 = [412,215,p_level1];
Plane1_P4 = [412,285,p_level1];
Plane1 = [Plane1_P1;Plane1_P2;Plane1_P3;Plane1_P4;Plane1_P1]; %last one is added for display
Plane1(:,2)=-Plane1(:,2); %flip due to image + is y-

Plane1_coor = zeros(5,3);
Plane1_coor(:,1) = (Plane1(:,1)-center_coor(1))*scale_bar;
Plane1_coor(:,2) = (Plane1(:,2)-center_coor(2))*scale_bar;
Plane1_coor(:,3) = (Plane1(:,3)-center_coor(3))*scale_bar;

Plane2_P1 = [517,216,p_level2];
Plane2_P2 = [517,244,p_level2];
Plane2_P3 = [517,287,p_level1];
Plane2_P4 = [517,216,p_level1];
Plane2 = [Plane2_P1;Plane2_P2;Plane2_P3;Plane2_P4;Plane2_P1];
Plane2(:,2)=-Plane2(:,2); %flip due to image + is y-

Plane2_coor = zeros(5,3);
Plane2_coor(:,1) = (Plane2(:,1)-center_coor(1))*scale_bar;
Plane2_coor(:,2) = (Plane2(:,2)-center_coor(2))*scale_bar;
Plane2_coor(:,3) = (Plane2(:,3)-center_coor(3))*scale_bar;

end

