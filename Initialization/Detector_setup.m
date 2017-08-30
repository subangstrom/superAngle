function [ Detector] = Detector_setup( Detector )
% Calculate x-ray with eular angles (in sphere coordinate) that enters
% the EDS detectors
% Detector parameters must be input from file
% Copyright by Weizong Xu, email to wxu4@ncsu.edu for any info and questions. 
% March 2015

tot_Det_num = Detector.tot_Det_num;
dphi=Detector.dAngle;
dtheta=Detector.dAngle;
%num_theta= double(int16(90.0/dAngle));
max_size=0;
angle_all=zeros(360/dphi*90/dtheta/4,5,tot_Det_num);

for i=1:tot_Det_num
    %Suggested detector parameters for input
    %Det_radius = 2.9; %Detector radius in mm (rd=2.9 from MM 20(2014)1318-1326))
    %Det_dist = 12; %Distance from Detector to eucentric point in mm (r), 15->solid angle for each detector is 0.13 if r=3.1
    %Det_takeoff = 18; %Take off angle in deg
    %Det_azimuth = 135, 45, 225, 315 deg; %%Azimuth angle in deg 90deg to the default thin wedge
    %Det_tilt = -7.5; %-(22.5-15)=-7.5; - more perpendicular to x-y %Detector tilt angle in deg
    %det_coor --> Center coordinate
    %D_n --> plane normal of the detector
    %e.g. [det_coor, D_n] = Get_Det_Coor(Det_dist, Det_takeoff, Det_azimuth, Det_tilt);

    %=======Reference========
    %   Geometry in Titan
    %           Y+
    %            |
    %       D1   |     D2  
    %     135deg |   45deg
    %           (Z)-----X+
    %            
    %       D4         D3
    %     225deg     315deg
    %            
    %========================
    %Det_label = detector_para.Dector(i);
    Det_dist = Detector.detector_para.Distance(i);
    Det_takeoff = Detector.detector_para.Takeoff_angle(i);
    Det_azimuth = Detector.detector_para.Azimuth_angle(i);
    Det_tilt = Detector.detector_para.Self_tilt(i);
    Det_radius = Detector.detector_para.Detector_radius(i);

    %*******Detector orientation calculation*******
    %[det_coor, D_n] = Get_Det_Coor(Det_dist, Det_takeoff, Det_azimuth, Det_tilt);
    theta = (90-Det_takeoff)*pi/180;
    phi = Det_azimuth*pi/180;
    dtilt = - Det_tilt * pi/180; %negative due to tilt angle difination
    det_coor= zeros(1,3);  %x-y-z coordinate
    det_coor(1) = Det_dist*sin(theta)*cos(phi);
    det_coor(2) = Det_dist*sin(theta)*sin(phi);
    det_coor(3) = Det_dist*cos(theta);
    
    %n=Det_dist*[1,theta+dtilt, phi];
    D_n= zeros(1,3);  %x-y-z coordinate
    D_n(1) = Det_dist*sin(theta+dtilt)*cos(phi);
    D_n(2) = Det_dist*sin(theta+dtilt)*sin(phi);
    D_n(3) = Det_dist*cos(theta+dtilt);
    D_n = D_n/norm(D_n);
    
    %**********Get Eular angle to detector*************
    %Angular search for Xray signal reach to one detector with parameter above
    %theta range from (0,pi/2)
    %phi range from (0,2pi)
    [Det_angle_search, det_disp, Detector_solid_angle] = Det_collect_search2(det_coor, D_n, dtheta, dphi, Det_radius);
    Detector.detector_para.Solid_angle(i)=Detector_solid_angle;
 
    temp_size = size(Det_angle_search);
    for m=1:temp_size(1)
        for n=1:8
            angle_all(m,n,i)=Det_angle_search(m,n);
        end
    end
    
    if (temp_size(1)>max_size)
        max_size = temp_size(1);
    end
    
    %angle_search(:,:,i)=Det_angle_search;
    
    if (i==1)
        det_disp_all=det_disp;
    else
        det_disp_all=[det_disp_all;det_disp];
    end

end    
    angle_search(:,:,:)=angle_all(1:max_size,:,:);
    Detector.angle_search=angle_search;

    % unpack --> angle_search_i = angle_search(:,:,i);

%end

    %Plot3Dim(det_disp_all, 'detector')
    
end

