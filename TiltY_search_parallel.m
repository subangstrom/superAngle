function [ Al_outY, Ni_outY] = TiltY_search_parallel( search_Deg, d_Deg, sample_para_in, holder_para, angle_search, Spurious)
%Calculate X-ray absorption ratio and absorption correction efficient
%Search Y tilt in range of -search_Deg to +search_Deg
% Weizong Xu, Feb. 2015

% Input: sample_para
% Predev_angle_X = sample_para_in(1); %deg Assume wedge line parallel to x-axis, thin area towards Y+
% Predev_angle_Y = sample_para_in(2);
% TiltX = sample_para_in(3);
% %TiltY = sample_para_in(4);
Thickness0 = sample_para_in.Thickness;
n= sample_para_in.Slice_t;
u1= sample_para_in.uA*1e-7;
u2= sample_para_in.uB*1e-7;
% HolderX = sample_para_in(9);
% HolderY = sample_para_in(10);
t_chk = sample_para_in.t_chk;
% cal_chk = sample_para_in(12);
% EleA_num = sample_para_in(13);
% EleA_shell = sample_para_in(14);
% EleB_num = sample_para_in(15);
% EleB_shell = sample_para_in(16);
% Atomic_ratio = sample_para_in(17);
% AB_Density = sample_para_in(18);
% k_factors_other = sample_para_in(19);
% k_AB_ideal = sample_para_in(20);
% DepthZ = sample_para_in(21);
% probe_Ne = sample_para_in(22);
% acquire_time = sample_para_in(23);
%sample_para=sample_para_in; %cannot use this in parfor

%dt=Thickness/n; % n Slices along thickness direction

tot_search= double(int16(search_Deg*2/d_Deg+1));

Al_outY=zeros(tot_search,2);
Ni_outY=zeros(tot_search,2);
%Absrp_outY=zeros(tot_search,2);
%Xray_num=zeros(tot_search,2);
p=ProgressBar(tot_search);
parfor i=1:tot_search  %parfor- parallel for loop

        TiltY=(i-1)*d_Deg-search_Deg;
        if (t_chk == 1)
            Thickness = Thickness0;
        else
            Thickness = Thickness0/cos(TiltX*pi/180)/cos(TiltY*pi/180);
        end
        dt=Thickness/n; 
        sample_para=sample_para_in;
        sample_para.TiltY=TiltY;        
%         sample_para =[Predev_angle_X, Predev_angle_Y, TiltX, TiltY, Thickness, n, u1, u2, HolderX, HolderY, t_chk, cal_chk, EleA_num, EleA_shell, EleB_num, EleB_shell, ...
%             Atomic_ratio, AB_Density, k_factors_other, k_AB_ideal, DepthZ, probe_Ne, acquire_time];
        [H_angle_out] = Holder_shadow(sample_para, holder_para, angle_search);
        temp_num=size(H_angle_out);              
      if (temp_num(1)>=1)
            
        S_n0 = [0,0,1];
        [S_n] = sample_normal(S_n0,sample_para);
        
        angle_precal=zeros(temp_num(1),3);
        angle_precal(:,1)=S_n(1).*sin(H_angle_out(:,1)).*cos(H_angle_out(:,2))+ ...
            S_n(2).*sin(H_angle_out(:,1)).*sin(H_angle_out(:,2))+S_n(3).*cos(H_angle_out(:,1));
        angle_precal(:,2)=H_angle_out(:,3); %for Al
        angle_precal(:,3)=H_angle_out(:,4); %for Ni

        [ tot_A1, tot_A2] = absorb_cal_ultrafast (u1,u2, angle_precal, dt, Thickness, S_n);
        
%        if (temp_num(1)>=1)        
%        Al_outY(i,1)=i; %must be avoid when use parfor
        Al_outY(i,2)=tot_A1/n;%+Spurious(1);  %avg Absorption Al

        Ni_outY(i,2)=tot_A2/n;%+Spurious(2);  %avg Absorption Ni

        else
            Al_outY(i,2)=0;
            Ni_outY(i,2)=1e-20;
       end
p.progress;
end
Al_outY(:,1)=-search_Deg:d_Deg:search_Deg;
Ni_outY(:,1)=-search_Deg:d_Deg:search_Deg;
p.stop;

end