function [ holder_para_out ] = update_holder_para( holder_para, sample_para )
%Design for GUI interface
%Weizong Xu, July 2015
    
Depth = holder_para(1);
Effective_ring_diameter = holder_para(2);
Cal_side_effect = holder_para(3);
Wall_angle = holder_para(4);
%uA_holder = holder_para(5);
%uB_holder = holder_para(6);
grid_chk = holder_para(7);
type_grid = holder_para(8);
open_diameter_grid = holder_para(9);
chk_holder = holder_para(10);
Be_chk = holder_para(11);

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

holder_para_out = [Depth, Effective_ring_diameter, Cal_side_effect, Wall_angle, uA_holder, uB_holder, grid_chk, type_grid, open_diameter_grid, chk_holder, Be_chk];



end

