function [ H_angle_out ] = Holder_carrier_shadow_up_2(sample_para, holder_para, angle_search_in)
%If holder carrier is made of Be, it will shadow low energy X-ray, but
%semi-transparent to high energy X-ray. This acts a kind of filtering.
%This code deals with such senario.
%Weizong Xu, April, 2015, wxu4@ncsu.edu

TiltX = sample_para(3)*pi/180;
TiltY = sample_para(4)*pi/180; 
HolderX = sample_para(9);
HolderY = sample_para(10);


%depth = sample_para(21)/(cos(TiltX)*cos(TiltY)); %e.g. correct depth increment during XY tilt
%H_coor = [0, 0, depth];

%H_center = [HolderX, HolderY, depth];
depth = sample_para(21);
H_size = holder_para(2)*0.5; %e.g. 2.5*0.5 mm
chk = holder_para(3);
Alpha = holder_para(4)*pi/180; %Wall_angle convert to degree
uAl_K_Be = holder_para(5);
uNi_K_Be = holder_para(6);

k=1/tan(Alpha);
Cone_P(1)=-HolderX;
Cone_P(2)=-HolderY;
Cone_P(3)=depth-tan(Alpha)*H_size;

%if (chk>0)


RotX= [1, 0, 0; 0, cos(TiltX), -sin(TiltX); 0, sin(TiltX), cos(TiltX)]; %Clockwise
RotY= [cos(TiltY), 0, sin(TiltY); 0, 1, 0; -sin(TiltY), 0, cos(TiltY)]; %Clockwise
RotX_rev= [1, 0, 0; 0, cos(-TiltX), -sin(-TiltX); 0, sin(-TiltX), cos(-TiltX)]; %Anti-Clockwise
RotY_rev= [cos(-TiltY), 0, sin(-TiltY); 0, 1, 0; -sin(-TiltY), 0, cos(-TiltY)]; %Anti-Clockwise
RotYX=RotY*RotX;
RotXY_rev=RotX_rev*RotY_rev;
H_n= [0, 0, 1]*RotYX;
S_n_R= [1,0,0]*RotYX;
S_n_L= [-1,0,0]*RotYX;
H_center = [-HolderX, -HolderY, depth]*RotYX;
%set up FEI holder, hex wall of holder carrier
carrier_width = 2.2; %distance from center to the side wall.
pre_coor_factor = carrier_width/sqrt(2);
C_wall_coor = zeros(1,3,8);
C_wall_coor(:,:,1) = [pre_coor_factor-HolderX,pre_coor_factor-HolderY, 0]*RotYX;
C_wall_coor(:,:,2) = [carrier_width-HolderX, 0-HolderY, 0]*RotYX;
C_wall_coor(:,:,3) = [pre_coor_factor-HolderX,-pre_coor_factor-HolderY, 0]*RotYX;
C_wall_coor(:,:,4) = [0-HolderX, -carrier_width-HolderY, 0]*RotYX;
C_wall_coor(:,:,5) = [-pre_coor_factor-HolderX,-pre_coor_factor-HolderY, 0]*RotYX;
C_wall_coor(:,:,6) = [-carrier_width-HolderX, 0-HolderY, 0]*RotYX;
C_wall_coor(:,:,7) = [-pre_coor_factor-HolderX,pre_coor_factor-HolderY, 0]*RotYX;
C_wall_coor(:,:,8) = [0-HolderX, carrier_width-HolderY, 0]*RotYX;
C_wall_n = C_wall_coor;

%Initial coordinate of clip side legs to calculate clip blocking
d_cal=-0.05; %clip leg position relative to carrier plane (in mm)
depthZ=depth+d_cal; % size of the clip along z direction relative to the carrier front surface.
%along x-y plane  y=-x+sqrt(2), P_middle=carrier_width/sqrt(2),+-halfwidth/sqrt(2)
Clip_width = 0.9; %in mm;
C_x1 = carrier_width/sqrt(2)+Clip_width/2/sqrt(2);
C_y1 = carrier_width/sqrt(2)-Clip_width/2/sqrt(2);
C_x2 = C_y1;
C_y2 = C_x1;

P1a = [C_x1-HolderX,C_y1-HolderY,depthZ]*RotYX; %@C_wall_coor(:,:,1) position
P1b = [C_x2-HolderX,C_y2-HolderY,depthZ]*RotYX; %Pa---Pb
P1c = [C_x2-HolderX,C_y2-HolderY,0]*RotYX;     %|     |
P1d = [C_x1-HolderX,C_y1-HolderY,0]*RotYX;     %Pd---Pc

P3a = [C_x2-HolderX,-C_y2-HolderY,depthZ]*RotYX;%@C_wall_coor(:,:,38) position
P3b = [C_x1-HolderX,-C_y1-HolderY,depthZ]*RotYX;
P3c = [C_x1-HolderX,-C_y1-HolderY,0]*RotYX;
P3d = [C_x2-HolderX,-C_y2-HolderY,0]*RotYX;

P1b_a = (P1b-P1a)';
P1d_a = (P1d-P1a)'; 
P1a_b = -P1b_a; %(P1a-P1b)'; 
P1c_b = (P1c-P1b)'; 
P1b_c = -P1c_b; %(P1b-P1c)'; 
P1d_c = (P1d-P1c)'; 
P1a_d = -P1d_a; %(P1a-P1d)'; 
P1c_d = -P1d_c; %(P1c-P1d)';

P3b_a = (P3b-P3a)';
P3d_a = (P3d-P3a)'; 
P3a_b = -P3b_a; %(P3a-P3b)'; 
P3c_b = (P3c-P3b)'; 
P3b_c = -P3c_b; %(P3b-P3c)'; 
P3d_c = (P3d-P3c)'; 
P3a_d = -P3d_a; %(P3a-P3d)'; 
P3c_d = -P3d_c; %(P3c-P3d)';

[ S_L, S_R ] = holder_carrier_FEI_up();
S_La = [S_L(1,1)-HolderX,S_L(1,2)-HolderY,S_L(1,3)+depthZ]*RotYX; %@Side_wall_on the left
S_L_shift(:,1)=S_L(:,1)-HolderX;
S_L_shift(:,2)=S_L(:,2)-HolderY;
S_L_shift(:,3)=S_L(:,3)+depthZ;
% S_Lb = [S_L(2,1)-HolderX,S_L(2,2)-HolderY,S_L(2,3)+depthZ]*RotYX; %Pa---Pb
% S_Lc = [S_L(3,1)-HolderX,S_L(3,2)-HolderY,S_L(3,3)+depthZ]*RotYX; %|     |
% S_Ld = [S_L(4,1)-HolderX,S_L(4,2)-HolderY,S_L(4,3)+depthZ]*RotYX; %Pd---Pc

S_Ra = [S_R(1,1)-HolderX,S_R(1,2)-HolderY,S_R(1,3)+depthZ]*RotYX;%@Side_wall_on the right
S_R_shift(:,1)=S_R(:,1)-HolderX;
S_R_shift(:,2)=S_R(:,2)-HolderY;
S_R_shift(:,3)=S_R(:,3)+depthZ;

temp=size(angle_search_in);
Angle_size=temp(1);

H_angle = zeros(Angle_size,4);
H_xyz = zeros(Angle_size,3);
AS_xyz = zeros(Angle_size,3);
n=0;
for i=1:Angle_size
    X_theta = angle_search_in(i,1);
    X_phi= angle_search_in(i,2);
    X_int_Al= angle_search_in(i,3);
    X_int_Ni= angle_search_in(i,4);
    Ang_x = sin(X_theta)*cos(X_phi);
    Ang_y = sin(X_theta)*sin(X_phi);
    Ang_z = cos(X_theta); 
    Xray_dist = (H_n(1)*H_center(1)+H_n(2)*H_center(2)+H_n(3)*H_center(3))/ ...
                (H_n(1)*Ang_x+H_n(2)*Ang_y+H_n(3)*Ang_z);
            
    %Xray_dist = (H_n(1)*H_coor(1)+H_n(2)*H_coor(2)+H_n(3)*H_coor(3))/ ...
    %            (H_n(1)*sin(X_theta)*cos(X_phi)+H_n(2)*sin(X_theta)*sin(X_phi)+H_n(3)*cos(X_theta));
    %similar to f(H_n,H_center)
%Consider X-ray fully shadowed by holder with tilt angle  
                  %test only (non-shadow area)
                  AS_xyz(i,1) = Xray_dist*Ang_x;
                  AS_xyz(i,2) = Xray_dist*Ang_y;
                  AS_xyz(i,3) = Xray_dist*Ang_z;
                  %*********
            if (Xray_dist>0)
                H_x = AS_xyz(i,1);
                H_y = AS_xyz(i,2);
                H_z = AS_xyz(i,3);
                holder_dist = sqrt((H_x-H_center(1))^2+(H_y-H_center(2))^2+(H_z-H_center(3))^2);
                if (holder_dist <= H_size)
                    n=n+1;
                    H_angle(n,1) = X_theta;
                    H_angle(n,2) = X_phi;
                    H_angle(n,3) = X_int_Al; %for Al;
                    H_angle(n,4) = X_int_Ni; %for Ni;
                    H_xyz(n,1) = H_x;
                    H_xyz(n,2) = H_y;
                    H_xyz(n,3) = H_z;
                end
                
                if (holder_dist > H_size && chk==2)

                    
                    %The following code is to calculate penetration depth 
                    %of an X-ray travel through the holder carrier.
                    %It is based on the analytical solution to get
                    %intercept point between cone surface and line
                    %The validation of this solution can also be tested in
                    %an approximate calculation list at the end of code section. 
                    %Cone shape is based on the wall angle of alpha and detph
                    %To simplify the calculation, lines are rotated back to zero tilt
                    
                    
                    % k=1/tan(Alpha); %They are moved to the front of the code
                    % Cone_P(1)=-HolderX;
                    % Cone_P(2)=-HolderY;
                    % Cone_P(3)=depth-tan(Alpha)*H_size;
                    %L_n(1)=H_x;
                    %L_n(2)=H_y;
                    %L_n(3)=H_z;

                    L_n = [H_x,H_y,H_z]*RotX_rev*RotY_rev;

                    A=L_n(1)^2+L_n(2)^2-k^2*L_n(3)^2;
                    B=2*L_n(1)*Cone_P(1)+2*L_n(2)*Cone_P(2)-2*k^2*L_n(3)*Cone_P(3);
                    C=Cone_P(1)^2+Cone_P(2)^2-k^2*Cone_P(3)^2;

                    %temp=B^2-4*A*C;
                    t1=(B+sqrt(B^2-4*A*C))/(2*A);
                    t2=(B-sqrt(B^2-4*A*C))/(2*A);

                    out_t1=[L_n(1)*t1,L_n(2)*t1,L_n(3)*t1];
                    out_t2=[L_n(1)*t2,L_n(2)*t2,L_n(3)*t2];

                    %if (temp<0)  %CHK
                    %B
                    %end

                    if (dot(L_n,out_t1) > 0 )
                    out_t = out_t1;
                    else
                        if (dot(L_n,out_t2) > 0 )
                        out_t = out_t2;
                        end
                    end
                    
%                     if ((B^2-4*A*C)<0 || (dot(L_n,out_t1) <= 0 && dot(L_n,out_t2) <= 0))
%                     dot(L_n,out_t1) 
%                     dot(L_n,out_t2)
%                     shadow_display_check(sample_para, holder_para, [X_theta,X_phi])
%                     end
                    holder_pene = Xray_dist-norm(out_t); 
                    %norm(out_t) is the distance from (000) to out_t(x,y,z),all light is set from (000)                    
                    
                    %The following calculation of penetration depth in the holder
                    %carrier is based on the angle relationship, which largely avoid
                    %solving equation to find the cross point of inner cone surface
                    %It is valid when HolderX=0 and HolderY=0, but it is an good approximation when 
                    %they are relatively small. 
                    
                    %holder_pene=(holder_dist-H_size)*(tan(Alpha)/sin(Beta))/(tan(Alpha)/tan(Beta)-1);
                    %Beta = pi/2-acos(dot([H_x,H_y,H_z],H_n)/(norm([H_x,H_y,H_z])*norm(H_n)));
                    %holder_pene_2=(holder_dist-H_size)*(cos(Beta)+sin(Beta)/tan(Alpha-Beta));

                    %validation
                    %if (abs(holder_pene_2-holder_pene)/abs(holder_pene)>0.1)
                    %    holder_pene
                    %    holder_pene_2
                    %end
                    
                    
                    Xray_dist_wall = 100; %just initital a big number
                    min_ii=0;
                    for ii=1:8
                        C_n = C_wall_n(:,:,ii);
                        C_coor = C_wall_coor(:,:,ii);
                        Xray_dist_temp = (C_n(1)*C_coor(1)+C_n(2)*C_coor(2)+C_n(3)*C_coor(3))/ ...
                        (C_n(1)*Ang_x+C_n(2)*Ang_y+C_n(3)*Ang_z);
                        %above, assume side wall is straight, calculate the distance from center to the wall
                        if (Xray_dist_temp > 0 && Xray_dist_temp < Xray_dist_wall)
                            Xray_dist_wall = Xray_dist_temp; %find the mininal distance 
                            min_ii=ii; % so that one of the 8 walls that this beam penetrate out can be known 
                        end
                    end
                    
                    %compare whether beam goes out from top surface or side walls 
                    %i.e. which distance (Xray_dist_wall and Xray_dist) is smallest
                    
                    
                   if (min_ii==1 || min_ii==8) %right wall
                       
                        Dist_Swall = (S_n_R(1)*S_Ra(1)+S_n_R(2)*S_Ra(2)+S_n_R(3)*S_Ra(3))/ ...
                        (S_n_R(1)*Ang_x+S_n_R(2)*Ang_y+S_n_R(3)*Ang_z);
           
                        S_wall_coor = Dist_Swall*[Ang_x,Ang_y,Ang_z];
                        S_wall_coor_N = S_wall_coor *RotXY_rev;
                        xvR= S_R_shift(:,2)'; %no x-coordinate, plane is perpendicular to x-axis
                        yvR= S_R_shift(:,3)';
                        xR=S_wall_coor_N(2);
                        yR=S_wall_coor_N(3);
                        chk_in=inpolygon(xR,yR,xvR,yvR);
                        if (chk_in==1) %Check if projection point is within the polygon
                            holder_pene = Xray_dist-Dist_Swall;
                        end
                   end
                    
                   if (min_ii==7 || min_ii==8) %left wall
                       
                        Dist_Swall = (S_n_L(1)*S_La(1)+S_n_L(2)*S_La(2)+S_n_L(3)*S_La(3))/ ...
                        (S_n_L(1)*Ang_x+S_n_L(2)*Ang_y+S_n_L(3)*Ang_z);
           
                        S_wall_coor = Dist_Swall*[Ang_x,Ang_y,Ang_z];
                        S_wall_coor_N = S_wall_coor *RotXY_rev;
                        xvL= S_L_shift(:,2)'; %no x-coordinate, plane is perpendicular to x-axis
                        yvL= S_L_shift(:,3)';
                        xL=S_wall_coor_N(2);
                        yL=S_wall_coor_N(3);
                        chk_in=inpolygon(xL,yL,xvL,yvL);
                        if (chk_in == 1) %Check if projection point is within the polygon
                            holder_pene = Xray_dist-Dist_Swall;
                        end
                   end
                    
                   
                    
                    if (Xray_dist < Xray_dist_wall) %beam goes out through top surface, not hit side wall
                        %note Xray_dist here is out of H_size of holes.

                    n=n+1;
                    H_angle(n,1) = X_theta;
                    H_angle(n,2) = X_phi;
                    H_angle(n,3) = X_int_Al*exp(-(uAl_K_Be*holder_pene/10)); %/10 since cm --> mm
                    H_angle(n,4) = X_int_Ni*exp(-(uNi_K_Be*holder_pene/10));
                    
                    else
                        
                    ii_chk=1;
                    Angle_chk=zeros(8,1);
                    %precal_ang=C_wall_n(1)*Ang_x+C_wall_n(2)*Ang_y+C_wall_n(3)*Ang_z;
                    if (min_ii==1)
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
%                        P_wall = Xray_dist_wall*[sin(X_theta)*cos(X_phi),sin(X_theta)*sin(X_phi),cos(X_theta)];
%                        Angle_chk(1)=acos(dot(P_wall-P1a,P1b-P1a)/(norm(P_wall-P1a)*norm(P1b-P1a)));
%                        Angle_chk(2)=acos(dot(P_wall-P1a,P1d-P1a)/(norm(P_wall-P1a)*norm(P1d-P1a))); 
%                        Angle_chk(3)=acos(dot(P_wall-P1b,P1a-P1b)/(norm(P_wall-P1b)*norm(P1a-P1b))); 
%                        Angle_chk(4)=acos(dot(P_wall-P1b,P1c-P1b)/(norm(P_wall-P1b)*norm(P1c-P1b))); 
%                        Angle_chk(5)=acos(dot(P_wall-P1c,P1b-P1c)/(norm(P_wall-P1c)*norm(P1b-P1c))); 
%                        Angle_chk(6)=acos(dot(P_wall-P1c,P1d-P1c)/(norm(P_wall-P1c)*norm(P1d-P1c))); 
%                        Angle_chk(7)=acos(dot(P_wall-P1d,P1a-P1d)/(norm(P_wall-P1d)*norm(P1a-P1d))); 
%                        Angle_chk(8)=acos(dot(P_wall-P1d,P1c-P1d)/(norm(P_wall-P1d)*norm(P1c-P1d)));
%                        if (max(Angle_chk)<pi/2) %Check if projection point is within Pa-Pd
                            ii_chk=0; %beam is blocked by the clip at wall #1
                        %for test purpose
                        %K=[P1a;P1b;P1c;P1d];
                        %scatter3(K(:,1),K(:,2),K(:,3));
                        %hold on;
                        %scatter3(P_wall(1),P_wall(2),P_wall(3),'r')
                        %hold off;
                        end
                    end
                                        
                    if (min_ii==3)
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
%                        P_wall = Xray_dist_wall*[sin(X_theta)*cos(X_phi),sin(X_theta)*sin(X_phi),cos(X_theta)];
%                        Angle_chk(1)=acos(dot(P_wall-P3a,P3b-P3a)/(norm(P_wall-P3a)*norm(P3b-P3a)));
%                        Angle_chk(2)=acos(dot(P_wall-P3a,P3d-P3a)/(norm(P_wall-P3a)*norm(P3d-P3a))); 
%                        Angle_chk(3)=acos(dot(P_wall-P3b,P3a-P3b)/(norm(P_wall-P3b)*norm(P3a-P3b))); 
%                        Angle_chk(4)=acos(dot(P_wall-P3b,P3c-P3b)/(norm(P_wall-P3b)*norm(P3c-P3b))); 
%                        Angle_chk(5)=acos(dot(P_wall-P3c,P3b-P3c)/(norm(P_wall-P3c)*norm(P3b-P3c))); 
%                        Angle_chk(6)=acos(dot(P_wall-P3c,P3d-P3c)/(norm(P_wall-P3c)*norm(P3d-P3c))); 
%                        Angle_chk(7)=acos(dot(P_wall-P3d,P3a-P3d)/(norm(P_wall-P3d)*norm(P3a-P3d))); 
%                        Angle_chk(8)=acos(dot(P_wall-P3d,P3c-P3d)/(norm(P_wall-P3d)*norm(P3c-P3d)));
%                        if (max(Angle_chk)<pi/2) %Check if projection point is within Pa-Pd
                            ii_chk=0; %beam is blocked by the clip at wall #3
                        end
                    end
                    %if (phi_chk1>pi/20 && phi_chk2>pi/20) %simplified need construction here
                    %ii_chk=1;%disable clip wall blocking
                    if (ii_chk==1)
                    n=n+1;
                    H_angle(n,1) = X_theta;
                    H_angle(n,2) = X_phi;
                    
                    
                    C_x = Xray_dist_wall*Ang_x;
                    C_y = Xray_dist_wall*Ang_y;
                    C_z = Xray_dist_wall*Ang_z;
                    wall_pene = holder_pene-sqrt((C_x-H_x)^2+(C_y-H_y)^2+(C_z-H_z)^2);
                    H_angle(n,3) = X_int_Al*exp(-(uAl_K_Be*wall_pene/10)); %/10 since cm --> mm
                    H_angle(n,4) = X_int_Ni*exp(-(uNi_K_Be*wall_pene/10));                       
                    
                    end
                    %H_angle(n,3) = X_int_Al*exp(-(uAl_K_Be*holder_pene/10)); %/10 since cm --> mm
                    %H_angle(n,4) = X_int_Ni*exp(-(uNi_K_Be*holder_pene/10));
                    
                    
                    %H_angle(n,3) = 0;
                    %H_angle(n,4) = 0;
                    %else
                    %H_angle(n,4) = X_int_Ni*exp(-(uNi_K_Be*wall_pene/10));
                    %end
                    %H_angle(n,5) = holder_pene;
                    end
                
                
                end %%end if chk==2
            end
end
H_angle_out = zeros(n,4);
data_disp = zeros(n,3);
for i=1:n
    H_angle_out(i,1)=H_angle(i,1);
    H_angle_out(i,2)=H_angle(i,2);
    H_angle_out(i,3)=H_angle(i,3);
    H_angle_out(i,4)=H_angle(i,4);
    %H_angle_out(i,5)=H_angle(i,5);
    data_disp(i,1)=H_xyz(i,1);
    data_disp(i,2)=H_xyz(i,2);
    data_disp(i,3)=H_xyz(i,3);
end
%Display for testing
%sample_para(3)
%Scatter3Dim(data_disp, 'holder shadow');
%n
%Scatter3Dim(AS_xyz, 'non-shadow');
%temp(1)
%sample_para

%else
%    H_angle_out=angle_search_in;
%end

end

