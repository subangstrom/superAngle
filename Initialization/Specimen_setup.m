function [ sample_para ] = Specimen_setup( Specimen_file, t_chk )
%Setup specimen parameters from file
%Weizong Xu, Feb. 2015, wxu4@ncsu.edu

%**********Setup specimen************

if exist(Specimen_file, 'file')
    Datainput = xlsread(Specimen_file);
else
    uiwait(msgbox('Error!! Unable to find input file for specimen'));
    uiwait(msgbox('Demo, Ni3Al sample with t=50nm, inclination_angle 0deg, calculate Al-K/Ni-K'));
    Datainput = [3,0,0,50,100,0,0,0,0,13,1,28,1,0,1/3,7.5,1,0.22,10e8,100];
end


sample_para.Deviation_angle_X = Datainput(1,2); %angle between top surface and x-axis of the holder
sample_para.Deviation_angle_Y = Datainput(1,3); %angle between top surface and y-axis of the holder
sample_para.Thickness = Datainput(1,4); %sample thickness in nm 
sample_para.Slice_t = Datainput(1,5); % Number of Slice along thickness direction for integration

%Initialize sample with zero tilt
sample_para.TiltX = Datainput(1,6); %Xtilt in deg, useful for Point_search
sample_para.TiltY = Datainput(1,7); %Ytilt in deg, useful for Point_search & TiltX_search
sample_para.POSX = Datainput(1,8); %mm
sample_para.POSY = Datainput(1,9);

%Wedge_angle = 10;
%TiltY = 0;
%RotateZ = 0;
%POSX =0;
%POSY= 0;

sample_para.t_chk=t_chk;

sample_para.EleA_num = Datainput(1,10);
sample_para.EleA_shell = Datainput(1,11); %1) K; 2) L; 3) M;
sample_para.EleB_num = Datainput(1,12);
sample_para.EleB_shell = Datainput(1,13); %1) K; 2) L; 3) M;
sample_para.cal_chk = Datainput(1,14); % 1) cal composition from exp data, others compare with ideal composition
sample_para.Atomic_ratio = Datainput(1,15);
sample_para.AB_Density= Datainput(1,16); %g/cm3
sample_para.k_factors_other= Datainput(1,17); %other factors (unknown, calibrated from exp) that changes k-ratio, 1=ideal
sample_para.DepthZ = Datainput(1,18); %ideal 0.22mm
sample_para.probe_Ne = Datainput(1,19);
sample_para.acquire_time = Datainput(1,20);

[uA, uB, k_AB_ideal] = read_element( sample_para );

sample_para.uA=uA;
sample_para.uB=uB;
sample_para.k_AB_ideal=k_AB_ideal;
sample_para.Z_rotation=0; %degree
sample_para.chk_2nd=-1; %1 consider 2nd interception (25-100x slow), other not consider

[ convert_factor_A,convert_factor_B, ~,~ ] = absolute_scale_factor( sample_para );
sample_para.convert_factor_A=convert_factor_A;
sample_para.convert_factor_B=convert_factor_B;
sample_para.p_tetra_num=[];
%absorption coefficient database for Ni3Al
%uAl_K = Datainput(1,10);%4311.95*7.26; %(u/p)Al-K_Ni3Al*p
%uNi_K = 4494*7.2;
%uNi_K = Datainput(1,11);%60.012*7.26;  %(u/p)Ni-K_Ni3Al*p
% Note: I=I0*exp(-ul)=I0*exp[-(u/p)*pl]
% For Ni3Al system
% p=7.2 g/cm3
% (u/p)Al-K_Ni3Al = 4311.95 cm2/g (u/p)*p=31046 cm-1
% (u/p)Ni-K_Ni3Al = 61.012 cm2/g (u/p)*p=439.3 cm-1


%sample_para =[Deviation_angle_X, Deviation_angle_Y, TiltX, TiltY, Thickness, Slice_t, ...
%    uA, uB, POSX, POSY, t_chk, cal_chk, EleA_num, EleA_shell, EleB_num, EleB_shell, ...
%    Atomic_ratio, AB_Density, k_factors_other, k_AB_ideal, DepthZ, probe_Ne, acquire_time];



end

