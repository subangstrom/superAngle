function [S_n] = sample_normal(S_n0,sample_para)
%Determine sample normal direction
%Weizong Xu, Feb. 2015, wxu4@ncsu.edu

%anglesearch [theta,phi]

%Set specimen
Deviation_angle_X = sample_para.Deviation_angle_X; %angle between top surface and x-axis of the holder
Deviation_angle_Y = sample_para.Deviation_angle_Y; %angle between top surface and y-axis of the holder

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

TiltX = sample_para.TiltX; %deg temp one, will incoporate into loops
TiltY = sample_para.TiltY;


%probe_arriving --> [0, 0, t];
%S_coor = [0, 0, t];
%S_n0= [0, 0, 1];

Deviation_angle_X=Deviation_angle_X*pi/180;
Deviation_angle_Y=Deviation_angle_Y*pi/180;
TiltX=TiltX*pi/180;
TiltY=TiltY*pi/180;

RotX= [1, 0, 0; 0, cos(TiltX), -sin(TiltX); 0, sin(TiltX), cos(TiltX)]; %Clockwise
RotY= [cos(TiltY), 0, sin(TiltY); 0, 1, 0; -sin(TiltY), 0, cos(TiltY)]; %Clockwise

Rot_dev_angle_forY= [1, 0, 0; 0, cos(Deviation_angle_Y), -sin(Deviation_angle_Y); 0, sin(Deviation_angle_Y), cos(Deviation_angle_Y)]; %at Y angle: tilt along x-direction Clockwise
Rot_dev_angle_forX= [cos(Deviation_angle_X), 0, sin(Deviation_angle_X); 0, 1, 0; -sin(Deviation_angle_X), 0, cos(Deviation_angle_X)]; %at X angle: tilt along y-direction Clockwise
Rot_dev= Rot_dev_angle_forY * Rot_dev_angle_forX;
S_n= S_n0*Rot_dev*RotY*RotX;

end

