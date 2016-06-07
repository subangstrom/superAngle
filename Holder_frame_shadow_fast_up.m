function [ H_angle_out ] = Holder_frame_shadow_fast_up(sample_para, holder_para, angle_search_in)
%Shadowing effect from holder contour (FEI low background holder)
%Materials for holder contour is brass
%Weizong Xu, April, 2015, wxu4@ncsu.edu


HolderX = sample_para.POSX;
HolderY = sample_para.POSY;
%HolderX = HolderX+0;%temp set as zero, later consider HolderX, HolderY calebration with Stage display value
%HolderY = HolderY+0;
TiltX = sample_para.TiltX*pi/180;
TiltY = sample_para.TiltY*pi/180; %Note: holder frame won't tilt in Y
%direction, but z-height relative to the holder will change with HolderX
%(if it is not zero)
depthZ = sample_para.DepthZ;
depthZ = depthZ-HolderX*sin(TiltY); %compensate Z_shift during Y-tilt
HolderX = HolderX*cos(TiltY); %compensate X_shift during Y-tilt 
%Note: Need very accurate calibration of relative position between center point in STEM and holder
%center before use it. Also sample need stay at one point exactly. 
%If HolderX is small, simple results without compensation is also very close.
holder_frame_para=holder_para.FEI_frame_out;
FEI_frame_out = zeros(length(holder_frame_para),3);
FEI_frame_out(:,1) = holder_frame_para(:,1)-HolderX;
FEI_frame_out(:,2) = holder_frame_para(:,2)-HolderY;
FEI_frame_out(:,3) = holder_frame_para(:,3)+ depthZ;

FEI_frame_out_n(:,1) = holder_frame_para(:,4);
FEI_frame_out_n(:,2) = holder_frame_para(:,5);
FEI_frame_out_n(:,3) = holder_frame_para(:,6);
%distance_out= zeros(length(FEI_frame_out),1);
%distance_out(:,1) = sqrt((FEI_frame_out(:,1)).^2+(FEI_frame_out(:,2)).^2+(FEI_frame_out(:,3)).^2);


%Note: holder frame won't tilt in Y
RotX= [1, 0, 0; 0, cos(TiltX), -sin(TiltX); 0, sin(TiltX), cos(TiltX)]; %Clockwise
%HF_n_H= [0, 0, 1]*RotX;
%HF_n_V= [1, 0, 0]*RotX;
temp=size(angle_search_in);
Angle_size=temp(1);

Holder_frame = FEI_frame_out*RotX;
Holder_frame_normal = FEI_frame_out_n*RotX; 
HF_n_V=[1,0,0]*RotX;
C_wall_n=HF_n_V;



P1a = [holder_frame_para(1,8)-HolderX,holder_frame_para(2,8)-HolderY,holder_frame_para(3,8)+depthZ]*RotX;   %@holder wall at D1 path
P1b = [holder_frame_para(4,8)-HolderX,holder_frame_para(5,8)-HolderY,holder_frame_para(6,8)+depthZ]*RotX;   %Pa---Pb
P1c = [holder_frame_para(7,8)-HolderX,holder_frame_para(8,8)-HolderY,holder_frame_para(9,8)+depthZ]*RotX;   %|     |
P1d = [holder_frame_para(10,8)-HolderX,holder_frame_para(11,8)-HolderY,holder_frame_para(12,8)+depthZ]*RotX;%Pd---Pc

P2a = [holder_frame_para(13,8)-HolderX,holder_frame_para(14,8)-HolderY,holder_frame_para(15,8)+depthZ]*RotX;%@holder wall at D2 path
P2b = [holder_frame_para(16,8)-HolderX,holder_frame_para(17,8)-HolderY,holder_frame_para(18,8)+depthZ]*RotX;%Pa---Pb
P2c = [holder_frame_para(19,8)-HolderX,holder_frame_para(20,8)-HolderY,holder_frame_para(21,8)+depthZ]*RotX;%|     |
P2d = [holder_frame_para(22,8)-HolderX,holder_frame_para(23,8)-HolderY,holder_frame_para(24,8)+depthZ]*RotX;%Pd---Pc

P3a = [holder_frame_para(25,8)-HolderX,holder_frame_para(26,8)-HolderY,holder_frame_para(27,8)+depthZ]*RotX;%@holder wall at D3 path
P3b = [holder_frame_para(28,8)-HolderX,holder_frame_para(29,8)-HolderY,holder_frame_para(30,8)+depthZ]*RotX;%Pa---Pb
P3c = [holder_frame_para(31,8)-HolderX,holder_frame_para(32,8)-HolderY,holder_frame_para(33,8)+depthZ]*RotX;%|     |
P3d = [holder_frame_para(34,8)-HolderX,holder_frame_para(35,8)-HolderY,holder_frame_para(36,8)+depthZ]*RotX;%Pd---Pc

P4a = [holder_frame_para(37,8)-HolderX,holder_frame_para(38,8)-HolderY,holder_frame_para(39,8)+depthZ]*RotX;%@holder wall at D4 path
P4b = [holder_frame_para(40,8)-HolderX,holder_frame_para(41,8)-HolderY,holder_frame_para(42,8)+depthZ]*RotX;%Pa---Pb
P4c = [holder_frame_para(43,8)-HolderX,holder_frame_para(44,8)-HolderY,holder_frame_para(45,8)+depthZ]*RotX;%|     |
P4d = [holder_frame_para(46,8)-HolderX,holder_frame_para(47,8)-HolderY,holder_frame_para(48,8)+depthZ]*RotX;%Pd---Pc

precal1=C_wall_n(1)*P1a(1)+C_wall_n(2)*P1a(2)+C_wall_n(3)*P1a(3);
precal2=C_wall_n(1)*P2a(1)+C_wall_n(2)*P2a(2)+C_wall_n(3)*P2a(3);
precal3=C_wall_n(1)*P3a(1)+C_wall_n(2)*P3a(2)+C_wall_n(3)*P3a(3);
precal4=C_wall_n(1)*P4a(1)+C_wall_n(2)*P4a(2)+C_wall_n(3)*P4a(3);

P1b_a = (P1b-P1a)';
P1d_a = (P1d-P1a)'; 
P1a_b = -P1b_a; %(P1a-P1b)'; 
P1c_b = (P1c-P1b)'; 
P1b_c = -P1c_b; %(P1b-P1c)'; 
P1d_c = (P1d-P1c)'; 
P1a_d = -P1d_a; %(P1a-P1d)'; 
P1c_d = -P1d_c; %(P1c-P1d)';

P2b_a = (P2b-P2a)';
P2d_a = (P2d-P2a)'; 
P2a_b = -P2b_a; %(P2a-P2b)'; 
P2c_b = (P2c-P2b)'; 
P2b_c = -P2c_b; %(P2b-P2c)'; 
P2d_c = (P2d-P2c)'; 
P2a_d = -P2d_a; %(P2a-P2d)'; 
P2c_d = -P2d_c; %(P2c-P2d)';

P3b_a = (P3b-P3a)';
P3d_a = (P3d-P3a)'; 
P3a_b = -P3b_a; %(P3a-P3b)'; 
P3c_b = (P3c-P3b)'; 
P3b_c = -P3c_b; %(P3b-P3c)'; 
P3d_c = (P3d-P3c)'; 
P3a_d = -P3d_a; %(P3a-P3d)'; 
P3c_d = -P3d_c; %(P3c-P3d)';

P4b_a = (P4b-P4a)';
P4d_a = (P4d-P4a)'; 
P4a_b = -P4b_a; %(P4a-P4b)'; 
P4c_b = (P4c-P4b)'; 
P4b_c = -P4c_b; %(P4b-P4c)'; 
P4d_c = (P4d-P4c)'; 
P4a_d = -P4d_a; %(P4a-P4d)'; 
P4c_d = -P4d_c; %(P4c-P4d)';

frame_size = length(Holder_frame);

H_angle = zeros(Angle_size,4);
n=0;
%****************Speed up pre-search****************
search_range = ceil(frame_size/6);
i_test=ceil(Angle_size/2);
if (i_test>=1)
    X_theta = angle_search_in(i_test,1);
    X_phi= angle_search_in(i_test,2);
    Ang_x = sin(X_theta)*cos(X_phi);
    Ang_y = sin(X_theta)*sin(X_phi);
    Ang_z = cos(X_theta);  
    
   angle_chk_max = -1;  %initial a small number, =acos(pi)  
   j_min=0; 
   for j=1:frame_size
        
        HF_x = Holder_frame(j,1);
        HF_y = Holder_frame(j,2);
        HF_z = Holder_frame(j,3);

        if (HF_x*Ang_x > 0 && HF_y*Ang_y > 0)
            angle_swap=(HF_x*Ang_x+HF_y*Ang_y+HF_z*Ang_z)/...
                sqrt((HF_x^2+HF_y^2+HF_z^2));
             if (angle_swap > angle_chk_max) %angle_swap is positive
                angle_chk_max = angle_swap;
                j_min = j;
            end
        end       
        
   end 
end
%*************end*********


for i=1:Angle_size
    X_theta = angle_search_in(i,1);
    X_phi= angle_search_in(i,2);
    X_int_Al= angle_search_in(i,3);
    X_int_Ni= angle_search_in(i,4);
    Ang_x = sin(X_theta)*cos(X_phi);
    Ang_y = sin(X_theta)*sin(X_phi);
    Ang_z = cos(X_theta);  
    
    %*********pre-screen the x-ray hit on the four side wall of holders*********
    precal_ang=C_wall_n(1)*Ang_x+C_wall_n(2)*Ang_y+C_wall_n(3)*Ang_z;
    chk_hit_wall=1;
    Angle_chk=zeros(8,1);

    Xray_dist_wall = precal1/(C_wall_n(1)*Ang_x+C_wall_n(2)*Ang_y+C_wall_n(3)*Ang_z);
    if (Xray_dist_wall>0)
        P_wall = Xray_dist_wall*[Ang_x,Ang_y,Ang_z];
        Angle_chk(1)=(P_wall-P1a)*P1b_a;
        Angle_chk(2)=(P_wall-P1a)*P1d_a; 
        Angle_chk(3)=(P_wall-P1b)*P1a_b; 
        Angle_chk(4)=(P_wall-P1b)*P1c_b; 
        Angle_chk(5)=(P_wall-P1c)*P1b_c; 
        Angle_chk(6)=(P_wall-P1c)*P1d_c; 
        Angle_chk(7)=(P_wall-P1d)*P1a_d; 
        Angle_chk(8)=(P_wall-P1d)*P1c_d;
        if (min(Angle_chk)>=0) %Check if projection point is within Pa-Pd
            chk_hit_wall=0; %beam is blocked by the clip at wall #1
        end
    end
    
    if (chk_hit_wall==1)
        Xray_dist_wall = precal2/precal_ang;
        if (Xray_dist_wall>0)
            P_wall = Xray_dist_wall*[Ang_x,Ang_y,Ang_z];
            Angle_chk(1)=(P_wall-P2a)*P2b_a;
            Angle_chk(2)=(P_wall-P2a)*P2d_a; 
            Angle_chk(3)=(P_wall-P2b)*P2a_b; 
            Angle_chk(4)=(P_wall-P2b)*P2c_b; 
            Angle_chk(5)=(P_wall-P2c)*P2b_c; 
            Angle_chk(6)=(P_wall-P2c)*P2d_c; 
            Angle_chk(7)=(P_wall-P2d)*P2a_d; 
            Angle_chk(8)=(P_wall-P2d)*P2c_d;
            if (min(Angle_chk)>=0) %Check if projection point is within Pa-Pd
                chk_hit_wall=0; %beam is blocked by the clip at wall #2
            end
        end
    end
    
    if (chk_hit_wall==1)
        Xray_dist_wall = precal3/precal_ang;
        if (Xray_dist_wall>0)
            P_wall = Xray_dist_wall*[Ang_x,Ang_y,Ang_z];
            Angle_chk(1)=(P_wall-P3a)*P3b_a;
            Angle_chk(2)=(P_wall-P3a)*P3d_a; 
            Angle_chk(3)=(P_wall-P3b)*P3a_b; 
            Angle_chk(4)=(P_wall-P3b)*P3c_b; 
            Angle_chk(5)=(P_wall-P3c)*P3b_c; 
            Angle_chk(6)=(P_wall-P3c)*P3d_c; 
            Angle_chk(7)=(P_wall-P3d)*P3a_d; 
            Angle_chk(8)=(P_wall-P3d)*P3c_d;
            if (min(Angle_chk)>=0) %Check if projection point is within Pa-Pd
                chk_hit_wall=0; %beam is blocked by the clip at wall #3
            end
        end
    end
    
    if (chk_hit_wall==1)
        Xray_dist_wall = precal4/precal_ang;
        if (Xray_dist_wall>0)
            P_wall = Xray_dist_wall*[Ang_x,Ang_y,Ang_z];
            Angle_chk(1)=(P_wall-P4a)*P4b_a;
            Angle_chk(2)=(P_wall-P4a)*P4d_a; 
            Angle_chk(3)=(P_wall-P4b)*P4a_b; 
            Angle_chk(4)=(P_wall-P4b)*P4c_b; 
            Angle_chk(5)=(P_wall-P4c)*P4b_c; 
            Angle_chk(6)=(P_wall-P4c)*P4d_c; 
            Angle_chk(7)=(P_wall-P4d)*P4a_d; 
            Angle_chk(8)=(P_wall-P4d)*P4c_d;
            if (min(Angle_chk)>=0) %Check if projection point is within Pa-Pd
                chk_hit_wall=0; %beam is blocked by the clip at wall #4
            end
        end
    end
%***************************************************************************

 if (chk_hit_wall==1)

    %angle_min=pi; %initial a large number
    %HF_dist= 100;
    HF_x_min = 0;
    HF_y_min = 0;
    HF_z_min = 0;
    HF_j_min = 0;

    angle_chk_max = -1;  %initial a small number, =acos(pi)  
    %for jj=1:frame_size
    for jj=j_min-search_range:1:j_min+search_range %shorten the search range
        j=jj;
        if (jj<1)
            j=jj+frame_size;
        end
        
        if (jj>frame_size)
            j=jj-frame_size;
        end
        
        HF_x = Holder_frame(j,1);
        HF_y = Holder_frame(j,2);
        HF_z = Holder_frame(j,3);

        if (HF_x*Ang_x > 0 && HF_y*Ang_y > 0)
           
            %angle_swap=acos((HF_x*Ang_x+HF_y*Ang_y+HF_z*Ang_z)/...
            %    sqrt((HF_x^2+HF_y^2+HF_z^2)*(Ang_x^2+Ang_y^2+Ang_z^2)));
            angle_swap=(HF_x*Ang_x+HF_y*Ang_y+HF_z*Ang_z)/...
                sqrt((HF_x^2+HF_y^2+HF_z^2));
                        %Ang_x^2+Ang_y^2+Ang_z^2 == 1
 %           if (angle_swap < angle_min)
             if (angle_swap > angle_chk_max) %angle_swap is positive
                angle_chk_max = angle_swap;
                %HF_dist = distance_out(j);
                HF_x_min = HF_x;
                HF_y_min = HF_y;
                HF_z_min = HF_z;
                HF_j_min = j;
            end
        end       
        
    end
    
    
    HF_n(1) = Holder_frame_normal(HF_j_min,1);
    HF_n(2) = Holder_frame_normal(HF_j_min,2);
    HF_n(3) = Holder_frame_normal(HF_j_min,3);
    Xray_dist = (HF_n(1)*HF_x_min+HF_n(2)*HF_y_min+HF_n(3)*HF_z_min)/ ...
                (HF_n(1)*sin(X_theta)*cos(X_phi)+HF_n(2)*sin(X_theta)*sin(X_phi)+HF_n(3)*cos(X_theta));   
    %P_x = Ang_x*Xray_dist;
    %P_y = Ang_y*Xray_dist;
    P_z = Ang_z*Xray_dist;
            
    if ( P_z>=Holder_frame(HF_j_min,3) );
        
          if (holder_frame_para(HF_j_min,7)==3)
            Xray_dist2 = (HF_n_V(1)*HF_x_min+HF_n_V(2)*HF_y_min+HF_n_V(3)*HF_z_min)/ ...
                (HF_n_V(1)*sin(X_theta)*cos(X_phi)+HF_n_V(2)*sin(X_theta)*sin(X_phi)+HF_n_V(3)*cos(X_theta));   
            P_x2 = Ang_x*Xray_dist2;
            %P_y2 = Ang_y*Xray_dist2;
            %P_z2 = Ang_z*Xray_dist2;
            
            if ( abs(P_x2)<Holder_frame(HF_j_min,1));
                n=n+1;
                H_angle(n,1) = X_theta;
                H_angle(n,2) = X_phi;
                H_angle(n,3) = X_int_Al; %for Al;
                H_angle(n,4) = X_int_Ni; %for Ni;
            end
          else
            n=n+1;
            H_angle(n,1) = X_theta;
            H_angle(n,2) = X_phi;
            H_angle(n,3) = X_int_Al; %for Al;
            H_angle(n,4) = X_int_Ni; %for Ni;

          end
    
%     figure(1);scatter3(Holder_frame(:,1),Holder_frame(:,2),Holder_frame(:,3))
%     hold on
%     scatter3(HF_x_min,HF_y_min,HF_z_min,200,'d','filled')
%     scatter3(0,0,0,200,'filled','d')
%     plot3([0,HF_x_min],[0,HF_y_min],[0,HF_z_min]);
%     plot3([0,Ang_x*3.5],[0,Ang_y*3.5,],[0,Ang_z*3.5]);
%     plot3([0,P_x],[0,P_y],[0,P_z]);
%     az = 30;
%     el = 18;
%     view(az, el);
%     hold off 

    else 
    
        if (holder_frame_para(HF_j_min,7)==2)
            Xray_dist = (HF_n_V(1)*HF_x_min+HF_n_V(2)*HF_y_min+HF_n_V(3)*HF_z_min)/ ...
                (HF_n_V(1)*sin(X_theta)*cos(X_phi)+HF_n_V(2)*sin(X_theta)*sin(X_phi)+HF_n_V(3)*cos(X_theta));   
            %P_x = Ang_x*Xray_dist;
            P_y = Ang_y*Xray_dist;
            %P_z = Ang_z*Xray_dist;
            
            if ( abs(P_y)<abs(Holder_frame(HF_j_min,2)));
                n=n+1;
                H_angle(n,1) = X_theta;
                H_angle(n,2) = X_phi;
                H_angle(n,3) = X_int_Al; %for Al;
                H_angle(n,4) = X_int_Ni; %for Ni;
           
%                 figure(3);scatter3(Holder_frame(:,1),Holder_frame(:,2),Holder_frame(:,3))
%                 hold on
%                 scatter3(HF_x_min,HF_y_min,HF_z_min,200,'d','filled')
%                 scatter3(0,0,0,200,'filled','d')
%                 plot3([0,HF_x_min],[0,HF_y_min],[0,HF_z_min]);
%                 plot3([0,Ang_x*3.5],[0,Ang_y*3.5,],[0,Ang_z*3.5]);
%                 plot3([0,P_x],[0,P_y],[0,P_z]);
%                 az = 30;
%                 el = 18;
%                 view(az, el);
%                 hold off 
            end
        
        else
%         figure(2);scatter3(Holder_frame(:,1),Holder_frame(:,2),Holder_frame(:,3))
%         hold on
%         scatter3(HF_x_min,HF_y_min,HF_z_min,200,'d','filled')
%         scatter3(0,0,0,200,'filled','d')
%         plot3([0,HF_x_min],[0,HF_y_min],[0,HF_z_min]);
%         plot3([0,Ang_x*3.5],[0,Ang_y*3.5,],[0,Ang_z*3.5]);
%         plot3([0,P_x],[0,P_y],[0,P_z]);
%         hold off
        end

    end
 end


%         figure(4);scatter3(Holder_frame(:,1),Holder_frame(:,2),Holder_frame(:,3))
%         hold on
%         %scatter3(HF_x_min,HF_y_min,HF_z_min,200,'d','filled')
%         
%         scatter3(0,0,0,200,'filled','o')
%         %plot3([0,HF_x_min],[0,HF_y_min],[0,HF_z_min]);
%         plot3([0,Ang_x*3.5],[0,Ang_y*3.5,],[0,Ang_z*3.5]);
%         %plot3([0,P_x],[0,P_y],[0,P_z]);
%         Plane1_coor =[P1a;P1b;P1c;P1d;P1a];
%         Plane2_coor =[P2a;P2b;P2c;P2d;P2a];
%         Plane3_coor =[P3a;P3b;P3c;P3d;P3a];
%         Plane4_coor =[P4a;P4b;P4c;P4d;P4a];
%         plot3(Plane1_coor(:,1),Plane1_coor(:,2),Plane1_coor(:,3))
%         plot3(Plane2_coor(:,1),Plane2_coor(:,2),Plane2_coor(:,3))
%         plot3(Plane3_coor(:,1),Plane3_coor(:,2),Plane3_coor(:,3))
%         plot3(Plane4_coor(:,1),Plane4_coor(:,2),Plane4_coor(:,3))
%         hold off

end


 H_angle_out = zeros(n,4);
 for i=1:n
     H_angle_out(i,1)=H_angle(i,1);
     H_angle_out(i,2)=H_angle(i,2);
     H_angle_out(i,3)=H_angle(i,3);
     H_angle_out(i,4)=H_angle(i,4);
 end

end

