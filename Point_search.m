function [ Al_out, Ni_out] = Point_search( sample_para, holder_para, holder_frame_para, angle_search, Spurious)
%Calculate X-ray absorption ratio and absorption correction efficient
%Search for one point
%Weizong Xu, Feb. 2015


%Wedge_angle = sample_para_in(1); %deg Assume wedge line parallel to x-axis, thin area towards Y+
%RotateZ = sample_para_in(2);
TiltX = sample_para(3);
TiltY = sample_para(4);

Thickness0 = sample_para(5);
n= sample_para(6);
t_chk = sample_para(11);

if (t_chk == 1)
    Thickness = Thickness0;
else
    Thickness = Thickness0/cos(TiltX*pi/180)/cos(TiltY*pi/180);
end
dt=Thickness/n; % n Slices along thickness direction

tot_Al=0;
Al_out=0;
tot_Ni=0;
Ni_out=0;

[H_angle_out] = Holder_shadow(sample_para, holder_para, holder_frame_para, angle_search);
temp_num=size(H_angle_out);              
[S_n] = sample_normal([0,0,1], sample_para);
        
angle_precal=zeros(temp_num(1),3);
angle_precal(:,1)=S_n(1).*sin(H_angle_out(:,1)).*cos(H_angle_out(:,2))+ ...
S_n(2).*sin(H_angle_out(:,1)).*sin(H_angle_out(:,2))+S_n(3).*cos(H_angle_out(:,1));
angle_precal(:,2)=H_angle_out(:,3); %for Al
angle_precal(:,3)=H_angle_out(:,4); %for Ni

for tt=dt:dt:Thickness
    [tot_Al_tt, tot_Ni_tt] = absorb_cal_fast(sample_para, angle_precal, tt, S_n);
    tot_Al=tot_Al+tot_Al_tt;
    tot_Ni=tot_Ni+tot_Ni_tt;
end

if (temp_num(1)>=1)
    Al_out=tot_Al/n+Spurious(1);  %avg Absorption Al
    Ni_out=tot_Ni/n+Spurious(2);  %avg Absorption Ni
else
    Al_out=0;
    Ni_out=1e-20;
end

end