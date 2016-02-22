function [ sample_para_out ] = update_specimen_para( sample_para )
%Design for GUI interface
%Weizong Xu, July 2015

Deviation_angle_X = sample_para(1);
Deviation_angle_Y = sample_para(2);
TiltX = sample_para(3);
TiltY = sample_para(4);
Thickness = sample_para(5);
Slice_t = sample_para(6);
%uA = sample_para(7);
%uB = sample_para(8);
POSX = sample_para(9);
POSY = sample_para(10);
t_chk = sample_para(11);
cal_chk = sample_para(12);
EleA_num = sample_para(13);
EleA_shell = sample_para(14);
EleB_num = sample_para(15);
EleB_shell = sample_para(16);
Atomic_ratio = sample_para(17);
AB_Density = sample_para(18);
k_factors_other = sample_para(19);
%k_AB_ideal = sample_para(20);
DepthZ = sample_para(21);
probe_Ne = sample_para(22);
acquire_time = sample_para(23);

[uA, uB, k_AB_ideal] = read_element( EleA_num, EleA_shell, EleB_num, EleB_shell, Atomic_ratio, AB_Density, k_factors_other, cal_chk);

%absorption coefficient database for Ni3Al
%uAl_K = Datainput(1,10);%4311.95*7.26; %(u/p)Al-K_Ni3Al*p
%uNi_K = 4494*7.2;
%uNi_K = Datainput(1,11);%60.012*7.26;  %(u/p)Ni-K_Ni3Al*p
% Note: I=I0*exp(-ul)=I0*exp[-(u/p)*pl]
% For Ni3Al system
% p=7.2 g/cm3
% (u/p)Al-K_Ni3Al = 4311.95 cm2/g (u/p)*p=31046 cm-1
% (u/p)Ni-K_Ni3Al = 61.012 cm2/g (u/p)*p=439.3 cm-1

sample_para_out =[Deviation_angle_X, Deviation_angle_Y, TiltX, TiltY, Thickness, Slice_t, ...
    uA, uB, POSX, POSY, t_chk, cal_chk, EleA_num, EleA_shell, EleB_num, EleB_shell, ...
    Atomic_ratio, AB_Density, k_factors_other, k_AB_ideal, DepthZ, probe_Ne, acquire_time];


end

