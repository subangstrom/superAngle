function [ model_sample, model_ether ] = tilt_model( model_sample, model_ether, sample_para )
%Weizong Xu, August, 2017

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
%TiltZ = sample_para.Z_rotation;
Deviation_angle_X = sample_para.Deviation_angle_X; %angle between top surface and x-axis of the holder
Deviation_angle_Y = sample_para.Deviation_angle_Y; %angle between top surface and y-axis of the holder
% if (Deviation_angle_X~=0 || Deviation_angle_Y~=0)
%     disp('WARNING: Deviation angle (X/Y) is not zero! When using model as imput, they should be normally set as ZERO!')
% end
Deviation_angle_X=Deviation_angle_X*pi/180;
Deviation_angle_Y=Deviation_angle_Y*pi/180;
TiltX=TiltX*pi/180;
TiltY=TiltY*pi/180;
%TiltZ=TiltZ*pi/180;

RotX= [1, 0, 0; 0, cos(TiltX), -sin(TiltX); 0, sin(TiltX), cos(TiltX)]; %Clockwise
RotY= [cos(TiltY), 0, sin(TiltY); 0, 1, 0; -sin(TiltY), 0, cos(TiltY)]; %Clockwise
%RotZ= [cos(TiltZ), -sin(TiltZ), 0; sin(TiltZ), cos(TiltZ),0; 0, 0, 1]; %Clockwise
Rot_dev_angle_forY= [1, 0, 0; 0, cos(Deviation_angle_Y), -sin(Deviation_angle_Y); 0, sin(Deviation_angle_Y), cos(Deviation_angle_Y)]; %at Y angle: tilt along x-direction Clockwise
Rot_dev_angle_forX= [cos(Deviation_angle_X), 0, sin(Deviation_angle_X); 0, 1, 0; -sin(Deviation_angle_X), 0, cos(Deviation_angle_X)]; %at X angle: tilt along y-direction Clockwise
Rot_dev= Rot_dev_angle_forY * Rot_dev_angle_forX;
%S_n= S_n0*Rot_dev*RotY*RotX;
Rot_M=Rot_dev*RotY*RotX; %require right multiple p_n*Rot_M p_n 1x3
%Rot_M=RotZ*Rot_dev*RotY*RotX; %require right multiple p_n*Rot_M p_n 1x3
%Note: for complicated sample, set RotZ, not set Rot_dev; for simple flat
%sample, set Rot_dev, not RotZ.


%model.p is 3x1 must be model.p'
model_sample.p=(model_sample.p'*Rot_M)';
model_ether.p=(model_ether.p'*Rot_M)';

%look at it
% fmodel_sample.p=(model_sample.p'*Rot_M)';
% fmodel_ether.p=(model_ether.p'*Rot_M)';
% fmodel_sample.t=model_sample.t;
% fmodel_ether.t=model_ether.t;
% figure;
% % pdeplot3D(model_sample.p,model_sample.t)
% % hold on;
% % pdeplot3D(model_ether.p,model_ether.t)
% % hold on;
% pdeplot3D(fmodel_sample.p,fmodel_sample.t)
% hold on;
% pdeplot3D(fmodel_ether.p,fmodel_ether.t,'FaceAlpha',0)
% hold off;


end

