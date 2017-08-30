function [ model_sample, model_ether ] = tilt_model_precal( TiltX, TiltY, TiltZ, model_sample, model_ether)
%Systematically tilt of specimen during setup
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


TiltX=TiltX*pi/180;
TiltY=TiltY*pi/180;
TiltZ=TiltZ*pi/180;

RotX= [1, 0, 0; 0, cos(TiltX), -sin(TiltX); 0, sin(TiltX), cos(TiltX)]; %Clockwise
RotY= [cos(TiltY), 0, sin(TiltY); 0, 1, 0; -sin(TiltY), 0, cos(TiltY)]; %Clockwise
RotZ= [cos(TiltZ), -sin(TiltZ), 0; sin(TiltZ), cos(TiltZ),0; 0, 0, 1]; %Clockwise
Rot_M=RotZ*RotY*RotX; %require right multiple p_n*Rot_M p_n 1x3


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

