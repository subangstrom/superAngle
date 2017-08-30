function [ det_angle_search, det_disp, Detector_solid_angle ] = Det_collect_search2(det_coor, D_n, dtheta_in, dphi_in, Det_radius)
%Weizong Xu, Feb. 2015
%Search the signal reach to one detector with parameter above
%theta range from (0?pi/2)
%phi range from (0,2pi)

%theta_start=0;
theta_end=0.5*pi;
dtheta=dtheta_in*pi/180;
%dtheta=0.5*pi/num_theta;
%[theta_search] = get_theta(num_theta);
theta_search=dtheta:dtheta:theta_end;
a_temp=length(theta_search);
%a_temp=num_theta;

%phi_start=0;
phi_end=2*pi;
dphi=dphi_in*pi/180;
phi_search=dphi:dphi:phi_end;
b_temp=length(phi_search);

det_disp=zeros(a_temp*b_temp,5);

L_D = [0, 0, 0];
temp_n=1;
for i=1:a_temp
    for j=1:b_temp
       X_theta = theta_search(i);
       X_phi = phi_search(j); 
       Xray_r = (D_n(1)*det_coor(1)+D_n(2)*det_coor(2)+D_n(3)*det_coor(3))/(D_n(1)*sin(X_theta)*cos(X_phi)+D_n(2)*sin(X_theta)*sin(X_phi)+D_n(3)*cos(X_theta));
       L_D(1) = Xray_r*sin(X_theta)*cos(X_phi);
       L_D(2) = Xray_r*sin(X_theta)*sin(X_phi);
       L_D(3) = Xray_r*cos(X_theta);
       %data (i,j,1)=L_D(1);
       %data (i,j,2)=L_D(2);
       %data (i,j,3)=L_D(3);
       chk_dist= sqrt((L_D(1)-det_coor(1))^2+(L_D(2)-det_coor(2))^2+(L_D(3)-det_coor(3))^2);
       if (chk_dist< Det_radius)
            det_disp(temp_n,1)=L_D(1);
            det_disp(temp_n,2)=L_D(2);
            det_disp(temp_n,3)=L_D(3);
            det_angle_search(temp_n,1)=X_theta;
            det_angle_search(temp_n,2)=X_phi;
            det_angle_search(temp_n,3)=sin(X_theta)*dtheta*dphi; %Set initial intensity (d_solid angle) for Al, an arbitrary value
            det_angle_search(temp_n,4)=det_angle_search(temp_n,3); % for Ni, it is equal right now, could be modified later.
            det_angle_search(temp_n,5)=det_angle_search(temp_n,3); %Set d_solid angle for this beam 
            det_angle_search(temp_n,6)=sin(X_theta)*cos(X_phi);% convert to Cartesian coordinate for 3D model usage
            det_angle_search(temp_n,7)=sin(X_theta)*sin(X_phi);
            det_angle_search(temp_n,8)=cos(X_theta);
            %output(temp_n,3)=1; %Set initial intensity for Al, an arbitrary value
            %output(temp_n,4)=1; % for Ni, it is equal right now, could be modified later.
            temp_n=temp_n+1;
            
       end
    end
end


%Display
%Detector_solid_angle=temp_n/(2*num_theta*b_temp(2))*4*pi
Detector_solid_angle=sum(det_angle_search(:,5));
%disp (['Single detector solid angle = ',num2str(sum(det_angle_search(:,5)))]); 
%Plot3Dim(det_disp, 'detector')
end
