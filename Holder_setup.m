function [ holder_para, FEI_frame_out ] = Holder_setup( chk_shadow, Holder_file, sample_para )
%setup holder parameter
%Weizong Xu, March, 2015, wxu4@ncsu.edu


Datainput = Excel_input(Holder_file);
if (Datainput==0)
    uiwait(msgbox('Unable to find input file for holder, set as ideal parameter for FEI low background Be holder'));
    holder_para = [0.22,3.1,2,20,0,0,0,0,0];
    return;
end

chk_holder = Datainput(1);
Depth=Datainput(2); %mm holder+grid+sample_incline
Effective_ring_diameter = Datainput(3); %mm

Cal_side_effect = chk_shadow; %Datainput(4);
%>0 and <>2) consider holder shadowing with fully cutoff, no Xray outside cutoff distance
%2) consider holder shadowing with cutoff allowing high energy Xray penetration
%<=0 value) not consider holder shadowing effect
Wall_angle = Datainput(4); %deg

%Check if it is FEI low background holder
Be_chk = Datainput (5);
grid_chk = Datainput (6);
type_grid = Datainput (7);
open_diameter_grid = Datainput (8);


EleA_num = sample_para(13);
EleA_shell = sample_para(14); %1) K; 2) L; 3) M;
EleB_num = sample_para(15);
EleB_shell = sample_para(16); %1) K; 2) L; 3) M;
Absorp_table (:,:,1) = xlsread ('Absorption coefficient_K.xlsx'); %K-shell cm2/g
Absorp_table (:,:,2) = xlsread ('Absorption coefficient_L.xlsx'); %L-shell cm2/g
Absorp_table (:,:,3) = xlsread ('Absorption coefficient_M.xlsx'); %M-shell cm2/g
if (Be_chk==1)
    Be_density = 1.848;%g/cm3
    uA_holder = Absorp_table(4,EleA_num,EleA_shell)*Be_density; %atomic number of Be is 4
    uB_holder = Absorp_table(4,EleB_num,EleB_shell)*Be_density;
else
    Mo_density = 10.2;
    uA_holder = Absorp_table(42,EleA_num,EleA_shell)*Mo_density; %assume it is Mo (atom#42)
    uB_holder = Absorp_table(42,EleB_num,EleB_shell)*Mo_density;


    %uAl_K_Mo= Datainput(5);%1956*10.2; %Absroption due to holder material (Mo) for Al and Ni
    %uNi_K_Mo= Datainput(6);%191.9*10.2;
end

holder_para = [Depth, Effective_ring_diameter, Cal_side_effect, Wall_angle, uA_holder, uB_holder, grid_chk, type_grid, open_diameter_grid, chk_holder, Be_chk];

[ FEI_frame_out] = holder_contour_FEI_up;

end

