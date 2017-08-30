function [ H_angle_out ] = Holder_grid_shadow_up(sample_para, holder_para, angle_search_in)
%Shadowing effect from supporting grid
%Also pre-shadow impossible beam to speed up calculation
%Material is from input, either Mo or Cu
%Weizong Xu, April, 2015,wxu4@ncsu.edu;

grid_align_X=0;
grid_align_Y=0;
H_depth = sample_para.DepthZ-holder_para.Depth; %ideal depth for Titan 0.22mm
if (holder_para.grid_chk==1) %chk if there is a grid on sample
    depth=H_depth;
    %depth =0.05;
    %H_size = holder_para(9)*0.5; %grid radius e.g. 2.5*0.5 mm
    
    if (holder_para.type_grid~=1)
        H_size = holder_para.open_diameter_grid;
    else
        H_size = 1.0*0.5; %0.5 shortest distance of the 1x2 slot grid
    end

    if (H_depth>=0.2 || H_depth<=0)
        uiwait(msgbox('WARNING!! Unreasonable input of grid depth, please check holder input file'));
    end

else %if not, set a large grid to pre-shadow impossible beam(speed up cal)
    depth = sample_para.DepthZ;
    H_size = 10; %a large number than holder size
end
    


TiltX = sample_para.TiltX*pi/180;
TiltY = sample_para.TiltY*pi/180;
HolderX = sample_para.POSX-grid_align_X;
HolderY = sample_para.POSY-grid_align_Y;


%depth = sample_para(21)/(cos(TiltX)*cos(TiltY)); %e.g. correct depth increment during XY tilt
%H_coor = [0, 0, depth];

%H_center = [HolderX, HolderY, depth];
%depth = holder_para(8); %50um sample_para(21);
%H_size = holder_para(9)*0.5; %grid radius e.g. 2.5*0.5 mm
%depth = 0.05; %sample_para(21)-0.22; %50um sample_para(21);
%H_size = 0.75; %grid radius e.g. 2.5*0.5 mm

RotX= [1, 0, 0; 0, cos(TiltX), -sin(TiltX); 0, sin(TiltX), cos(TiltX)]; %Clockwise
RotY= [cos(TiltY), 0, sin(TiltY); 0, 1, 0; -sin(TiltY), 0, cos(TiltY)]; %Clockwise
RotX_rev= [1, 0, 0; 0, cos(-TiltX), -sin(-TiltX); 0, sin(-TiltX), cos(-TiltX)]; %Anti-Clockwise
RotY_rev= [cos(-TiltY), 0, sin(-TiltY); 0, 1, 0; -sin(-TiltY), 0, cos(-TiltY)]; %Anti-Clockwise
RotXY_rev=RotX_rev*RotY_rev;
RotYX=RotY*RotX;
H_n= [0, 0, 1]*RotYX;
H_center = [-HolderX, -HolderY, depth]*RotYX;
H_center_rev = [-HolderX, -HolderY, depth];
%H_n= [0, 0, 1]*RotX;
%H_center = [HolderX, HolderY, depth]*RotX;
%H_coor= [0, 0, depth]*RotY*RotX;

temp=size(angle_search_in);
Angle_size=temp(1);

H_angle = zeros(Angle_size,8);
H_xyz = zeros(Angle_size,3);
AS_xyz = zeros(Angle_size,3);
n=0;
for i=1:Angle_size
    X_theta = angle_search_in(i,1);
    X_phi= angle_search_in(i,2);
    X_int_Al= angle_search_in(i,3);
    X_int_Ni= angle_search_in(i,4);
    tmp1=angle_search_in(i,5);
    tmp2=angle_search_in(i,6);
    tmp3=angle_search_in(i,7);
    tmp4=angle_search_in(i,8);   
    Xray_dist = (H_n(1)*H_center(1)+H_n(2)*H_center(2)+H_n(3)*H_center(3))/ ...
                (H_n(1)*sin(X_theta)*cos(X_phi)+H_n(2)*sin(X_theta)*sin(X_phi)+H_n(3)*cos(X_theta));
            
%Consider X-ray fully shadowed by holder with tilt angle  
                  %test only (non-shadow area)
                  AS_xyz(i,1) = Xray_dist*sin(X_theta)*cos(X_phi);
                  AS_xyz(i,2) = Xray_dist*sin(X_theta)*sin(X_phi);
                  AS_xyz(i,3) = Xray_dist*cos(X_theta);
                  %*********
               if (Xray_dist>0)
                H_x = Xray_dist*sin(X_theta)*cos(X_phi);
                H_y = Xray_dist*sin(X_theta)*sin(X_phi);
                H_z = Xray_dist*cos(X_theta);
                holder_dist = sqrt((H_x-H_center(1))^2+(H_y-H_center(2))^2+(H_z-H_center(3))^2);
                
                if (holder_dist <= H_size)
                    n=n+1;
                    H_angle(n,1) = X_theta;
                    H_angle(n,2) = X_phi;
                    H_angle(n,3) = X_int_Al; %for Al;
                    H_angle(n,4) = X_int_Ni; %for Ni;
                    H_angle(n,5) = tmp1;
                    H_angle(n,6) = tmp2;
                    H_angle(n,7) = tmp3;
                    H_angle(n,8) = tmp4;
                    H_xyz(n,1) = H_x;
                    H_xyz(n,2) = H_y;
                    H_xyz(n,3) = H_z;
                else
                    
                   if (holder_dist <= 1.0 && holder_para.type_grid==1 )  %2.0*0.5; max open diameter of 1x2 slot
                        
                      Grid_pos = [H_x, H_y, H_z] * RotXY_rev; 
                      chk_pos = Grid_pos-H_center_rev;
                      %H_center_rev = [HolderX, HolderY, depth];
%                        if (abs(chk_pos(3)) > 0.01)
%                         chk_pos(3)
%                        end
                        
                       if (abs(chk_pos(1))<=0.575) %1.15/2
                           if (abs(chk_pos(2))<0.5)
                              n=n+1;
                            H_angle(n,1) = X_theta;
                            H_angle(n,2) = X_phi;
                            H_angle(n,3) = X_int_Al; %for Al;
                            H_angle(n,4) = X_int_Ni; %for Ni;
                            H_angle(n,5) = tmp1;
                            H_angle(n,6) = tmp2;
                            H_angle(n,7) = tmp3;
                            H_angle(n,8) = tmp4;
                            H_xyz(n,1) = H_x;
                            H_xyz(n,2) = H_y;
                            H_xyz(n,3) = H_z;
                           end
                       end
                       
                       if (chk_pos(1)<-0.575) %1.15/2,left circle
                           chk_dist_cir=sqrt((chk_pos(1)+0.575)^2+chk_pos(2)^2);
                           if (chk_dist_cir<0.515)
                                n=n+1;
                                H_angle(n,1) = X_theta;
                                H_angle(n,2) = X_phi;
                                H_angle(n,3) = X_int_Al; %for Al;
                                H_angle(n,4) = X_int_Ni; %for Ni;
                                H_angle(n,5) = tmp1;
                                H_angle(n,6) = tmp2;
                                H_angle(n,7) = tmp3;
                                H_angle(n,8) = tmp4;
                                H_xyz(n,1) = H_x;
                                H_xyz(n,2) = H_y;
                                H_xyz(n,3) = H_z;
                           end
                       end
                       
                       if (chk_pos(1)>0.575) %1.15/2,left circle
                           chk_dist_cir=sqrt((chk_pos(1)-0.575)^2+chk_pos(2)^2);
                           if (chk_dist_cir<0.515)
                                n=n+1;
                                H_angle(n,1) = X_theta;
                                H_angle(n,2) = X_phi;
                                H_angle(n,3) = X_int_Al; %for Al;
                                H_angle(n,4) = X_int_Ni; %for Ni;
                                H_angle(n,5) = tmp1;
                                H_angle(n,6) = tmp2;
                                H_angle(n,7) = tmp3;
                                H_angle(n,8) = tmp4;
                                H_xyz(n,1) = H_x;
                                H_xyz(n,2) = H_y;
                                H_xyz(n,3) = H_z;
                           end
                       end                       
                       
                   end
                    
                    
                    
                end
                
                
               end
end
H_angle_out = zeros(n,8);
data_disp = zeros(n,3);
for i=1:n
    H_angle_out(i,1)=H_angle(i,1);
    H_angle_out(i,2)=H_angle(i,2);
    H_angle_out(i,3)=H_angle(i,3);
    H_angle_out(i,4)=H_angle(i,4);
    H_angle_out(i,5)=H_angle(i,5);
    H_angle_out(i,6)=H_angle(i,6);
    H_angle_out(i,7)=H_angle(i,7);
    H_angle_out(i,8)=H_angle(i,8);
    data_disp(i,1)=H_xyz(i,1);
    data_disp(i,2)=H_xyz(i,2);
    data_disp(i,3)=H_xyz(i,3);
end


%else
    %no grid is on the sample
%    H_angle_out=angle_search_in; 
%end

%Display for testing
%sample_para(3)
%Scatter3Dim(data_disp, 'holder shadow');
%n
%Scatter3Dim(AS_xyz, 'non-shadow');
%temp(1)
%sample_para

%shadow_display_call( sample_para, holder_para, H_angle_out )

end

