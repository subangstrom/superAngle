function [ uA, uB, k_AB_ideal] = read_element( EleA_num, EleA_shell, EleB_num, EleB_shell, Atomic_ratio, AB_Density, k_factors_other, cal_chk)
% Input parameters for mass absorption coefficient and k factor calculation
% Ternary compounds calculation is under development. However, well known SrTiO3 is incoporated in. 
% Weizong Xu, March, 2015, email to wxu4@ncsu.edu for any info and questions.  


%e.g. Al-13 Ni-28
%Ni3Al;
%atomic ratio of element A : B
%Atomic_ratio = 1/3;
%EleA_num = 13;
%EleA_shell = 1; %1) K; 2) L; 3) M;
%EleB_num = 28;
%EleB_shell = 1; %1) K; 2) L; 3) M;
%AB_Density= 7.5; %g/cm3

Absorp_table (:,:,1) = Excel_input ('Absorption coefficient_K.xlsx'); %K-shell cm2/g
Absorp_table (:,:,2) = Excel_input ('Absorption coefficient_L.xlsx'); %L-shell cm2/g
Absorp_table (:,:,3) = Excel_input ('Absorption coefficient_M.xlsx'); %M-shell cm2/g

temp = Excel_input ('Atomic_info.xlsx');
Xray_table (:,:,1) = temp(1:3,:); % K-shell xray energy, crossection, flurence yield 
Xray_table (:,:,2) = temp(4:6,:); % L-shell xray energy, crossection, flurence yield 
Xray_table (:,:,3) = temp(7:9,:); % M-shell xray energy, crossection, flurence yield 

Atomic_weight_A = get_element_weight(EleA_num);
Atomic_weight_B = get_element_weight(EleB_num);

if (cal_chk ~= 1)

%Get composition
Atomic_percent_A = Atomic_ratio/(1+Atomic_ratio);
Atomic_percent_B = 1 - Atomic_percent_A;
Weight_percent_A = Atomic_percent_A*Atomic_weight_A/(Atomic_percent_A*Atomic_weight_A+Atomic_percent_B*Atomic_weight_B);
Weight_percent_B = 1 - Weight_percent_A;

%Get absorption coefficient
if (Absorp_table(EleA_num,EleA_num,EleA_shell) < 0 || Absorp_table(EleB_num,EleA_num,EleA_shell) < 0 || ...
        Absorp_table(EleA_num,EleB_num,EleB_shell) < 0 || Absorp_table(EleB_num,EleB_num,EleB_shell) < 0)
    disp('WARNING! Mass absorption coefficient is missing or incorrect. Please check Absorption coefficient_K/L/M.xlsx!');
    uiwait(msgbox('WARNING! Mass absorption coefficient is missing or incorrect. Please check Absorption coefficient_K/L/M.xlsx!'));
end
uA = Absorp_table(EleA_num,EleA_num,EleA_shell)*Weight_percent_A + Absorp_table(EleB_num,EleA_num,EleA_shell)*Weight_percent_B;
uA = uA * AB_Density;
uB = Absorp_table(EleA_num,EleB_num,EleB_shell)*Weight_percent_A + Absorp_table(EleB_num,EleB_num,EleB_shell)*Weight_percent_B;
uB = uB * AB_Density;

else
    disp('WARNING! Input composition is set as "not known", absorption coefficient is adjust to 0');
    uiwait(msgbox('WARNING! Input composition is set as "not known", absorption coefficient is adjust to 0'));
    uA=0;
    uB=0;
end

%SrTiO3 three element, by pass calculation WARNING:return/comment when it is done
%A-Sr B-Ti
if (EleA_num==38 && EleB_num==22 && abs(AB_Density-4.91)<0.3)
Weight_percent_A = 87.62/(87.62+47.87+16*3); %Sr
Weight_percent_B = 47.87/(87.62+47.87+16*3); %Ti
Weight_percent_C = 48/(87.62+47.87+16*3);%O
uA = Absorp_table(38,38,EleA_shell)*Weight_percent_A + Absorp_table(22,38,EleA_shell)*Weight_percent_B+Absorp_table(8,38,EleA_shell)*Weight_percent_C;
uB = Absorp_table(38,22,EleB_shell)*Weight_percent_A + Absorp_table(22,22,EleB_shell)*Weight_percent_B+Absorp_table(8,8,EleB_shell)*Weight_percent_C;
%AB_Density=4.81 nature-5.13 synthetic;%
uA = uA * AB_Density;
uB = uB * AB_Density;
disp('WARNING! Absorption calculation for Sr Ti O ternary element in SrTiO3!');
uiwait(msgbox('WARNING! Absorption calculation for Sr Ti O ternary element in SrTiO3!'));
end
%Get detector efficiency;
EleA_energy = Xray_table(1,EleA_num,EleA_shell); %keV
EleB_energy = Xray_table(1,EleB_num,EleB_shell); %keV
[ Det_eff_A, Det_eff_B ] = get_Det_eff( EleA_energy, EleB_energy, 'SDD_windowless_efficiency.xlsx' );
%Det_eff_A=get_Det_eff(EleA_energy);
%Det_eff_B=get_Det_eff(EleB_energy);


%Calculate ideal k-ratio
%At composition condition percentage eleA/eleB=1
%k_AB use function k=k_others*(Q_B*w_B*A_A*e_B)*(Q_A*w_A*A_B*e_A)

if ( Xray_table(2,EleB_num,EleB_shell) <= 0 || Xray_table(2,EleA_num,EleA_shell) <= 0 )
    disp('WARNING! Cross-section infomation is missing or incorrect. Counts/ratio calculation will be unavailable or incorrect. Please input it in Atomic_info.xlsx!');
    uiwait(msgbox('WARNING! Cross-section infomation is missing or incorrect. Counts/ratio calculation will be unavailable or incorrect. Please input it in Atomic_info.xlsx!'));
end

if ( Xray_table(3,EleB_num,EleB_shell) <= 0 || Xray_table(3,EleA_num,EleA_shell) <= 0 )
    disp('WARNING! X-ray fluorescence yield infomation is missing or incorrect. Counts/ratio calculation will be unavailable or incorrect. Please input it in Atomic_info.xlsx!');
    uiwait(msgbox('WARNING! X-ray fluorescence yield infomation is missing or incorrect. Counts/ratio calculation will be unavailable or incorrect. Please input it in Atomic_info.xlsx!'));
end

k_AB_ideal = k_factors_other*(Xray_table(2,EleB_num,EleB_shell)*Xray_table(3,EleB_num,EleB_shell)*Atomic_weight_A*Det_eff_B)/ ...
    (Xray_table(2,EleA_num,EleA_shell)*Xray_table(3,EleA_num,EleA_shell)*Atomic_weight_B*Det_eff_A);


end