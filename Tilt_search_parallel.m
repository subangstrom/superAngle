function [ Al_out, Ni_out, Absrp_out, Xray_num_out] = Tilt_search_parallel(search_Deg, d_Deg, sample_para_in, holder_para, angle_search, Spurious)
%Calculate X-ray absorption ratio and absorption correction efficient
%Search both X and Y tilt
% Weizong Xu, March. 2015

% Predev_angle_X = sample_para(1); %deg Assume wedge line parallel to x-axis, thin area towards Y+
% Predev_angle_Y = sample_para(2);
% TiltX = 0; %no use, arbitrary number;sample_para(3);
% %TiltY = sample_para(4);
Thickness = sample_para_in.Thickness;
n= sample_para_in.Slice_t;
% u1= sample_para_in.uA*1e-7;
% u2= sample_para_in.uB*1e-7;
% HolderX = sample_para(9);
% HolderY = sample_para(10);
% t_chk = sample_para_in.t_chk;
% cal_chk = sample_para(12);
% EleA_num = sample_para(13);
% EleA_shell = sample_para(14);
% EleB_num = sample_para(15);
% EleB_shell = sample_para(16);
% Atomic_ratio = sample_para(17);
% AB_Density = sample_para(18);
% k_factors_other = sample_para(19);
% k_AB_ideal = sample_para(20);
% DepthZ = sample_para(21);
% probe_Ne = sample_para(22);
% acquire_time = sample_para(23);

chk = 1; %set as Xtilt
dt=Thickness/n; % n Slices along thickness direction

tot_search= double(int16(search_Deg*2/d_Deg+1));
Al_out=zeros(tot_search,tot_search);
Ni_out=zeros(tot_search,tot_search);

p=ProgressBar(tot_search);
parfor j=1:tot_search
    TiltY = search_Deg-(j-1)*d_Deg; %process percentage indication
    sample_para=sample_para_in;
    sample_para.TiltY=TiltY; 
%     sample_para =[Predev_angle_X, Predev_angle_Y, TiltX, TiltY, Thickness, n, u1, u2, HolderX, HolderY, t_chk, cal_chk, EleA_num, EleA_shell, EleB_num, EleB_shell, ...
%             Atomic_ratio, AB_Density, k_factors_other, k_AB_ideal, DepthZ, probe_Ne, acquire_time];
    [Al_outX, Ni_outX] = TiltX_search(search_Deg, d_Deg, sample_para, holder_para, angle_search, Spurious, chk);
    Al_out(:,j)=Al_outX(:,2);
    Ni_out(:,j)=Ni_outX(:,2);
    p.progress;

end
p.stop;

Al_out=Al_out';
Ni_out=Ni_out';


end