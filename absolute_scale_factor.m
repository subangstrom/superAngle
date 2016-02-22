function [ factor_A,factor_B, pre_A, pre_B ] = absolute_scale_factor( sample_para )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
temp = Excel_input ('Atomic_info.xlsx');
Xray_table (:,:,1) = temp(1:3,:); % K-shell xray energy, crossection, flurence yield 
Xray_table (:,:,2) = temp(4:6,:); % L-shell xray energy, crossection, flurence yield 
Xray_table (:,:,3) = temp(7:9,:); % M-shell xray energy, crossection, flurence yield 

%EleA_num=13; EleA_shell=1;EleB_num=28;EleB_shell=1;Atomic_ratio=1/3;
EleA_num=sample_para(13);
EleA_shell=sample_para(14);
EleB_num=sample_para(15);
EleB_shell=sample_para(16);
Atomic_ratio=sample_para(17);
t=sample_para(5)*1e-9; %thickness
p=sample_para(18)*1e3; %density
Ne = sample_para(22);
tau = sample_para(23);

if (Ne<=1)
    disp('Warning: No electron current input for counts calculation!')
    msgbox('Warning: No electron current input for counts calculation!')
    Ne=0.01;
end
if (tau<=1)
    disp('Warning: No dwell time input for counts calculation!')
    msgbox('Warning: No dwell time input for counts calculation!')
    tau=0.01;
end

Atomic_weight_A = get_element_weight(EleA_num);
Atomic_weight_B = get_element_weight(EleB_num);
EleA_energy = Xray_table(1,EleA_num,EleA_shell); %keV
EleB_energy = Xray_table(1,EleB_num,EleB_shell); %keV
[ Det_eff_A, Det_eff_B ] = get_Det_eff( EleA_energy, EleB_energy, 'SDD_windowless_efficiency.xlsx' );
Atomic_percent_A = Atomic_ratio/(1+Atomic_ratio);
Atomic_percent_B = 1 - Atomic_percent_A;
Weight_percent_A = Atomic_percent_A*Atomic_weight_A/(Atomic_percent_A*Atomic_weight_A+Atomic_percent_B*Atomic_weight_B);
Weight_percent_B = 1 - Weight_percent_A;



Nv=6.0221413e23; %mol-1
Q_a= Xray_table(2,EleA_num,EleA_shell) * 1e-28; %m2
w_a= Xray_table(3,EleA_num,EleA_shell);% percentage unitless
a_a=1;%unit less relative transition probability should =1
M_a=Atomic_weight_A*1e-3; %kg mol-1
C_a=Weight_percent_A;%unit less
%p=7.5e3; %kg m-3
%t=84e-9; %nm
%Ne= 1.7237e9*0.44; %s-1, number of electrons per sec for (0.127nA probe); 0.44 is calibrated from CCD
%Ip=1.60217657e-19; %electron charge C;
Ip=1;
%tau=135;%s live time
De=Ne*Ip*tau;
%omega_a=0.171; %sr sr/4pi unitless
eff_a=Det_eff_A; %unitless (percentage)

Q_b= Xray_table(2,EleB_num,EleB_shell) * 1e-28; %m2
w_b= Xray_table(3,EleB_num,EleB_shell);% percentage unitless
a_b=1;%unit less relative transition probability should =1
M_b=Atomic_weight_B*1e-3; %kg mol-1
C_b=Weight_percent_B;%unit less
%omega_b=0.171; %sr sr/4pi unitless
eff_b=Det_eff_B; %unitless (percentage)

%k_AB_ideal = (Q_b*w_b*M_a*eff_b)/(Q_a*w_a*M_b*eff_a);



%K=Nv*Q_a*w_a*a_a*C_a*p*t*De*omega_a*eff_a/(M_a*4*pi);
%I_a=Nv*Q_a*w_a*a_a*C_a*p*t*De*omega_a*eff_a/(M_a*4*pi);
%I_b=Nv*Q_b*w_b*a_b*C_b*p*t*De*omega_b*eff_b/(M_b*4*pi);
%k_AB_ideal_test=I_b/I_a*C_a/C_b;
%counts_ratio=I_a/I_b;
factor_A=Nv*Q_a*w_a*a_a*C_a*p*t*De*eff_a/(M_a*4*pi);
factor_B=Nv*Q_b*w_b*a_b*C_b*p*t*De*eff_b/(M_b*4*pi);

pre_A=Nv*Q_a*w_a*a_a*p*t*De*eff_a/(M_a*4*pi);
pre_B=Nv*Q_b*w_b*a_b*p*t*De*eff_b/(M_b*4*pi);
end

