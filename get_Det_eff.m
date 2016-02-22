function [ eff_A, eff_B ] = get_Det_eff( energy_A, energy_B, filename )
%Get detector efficiency for Bruker Windowless SDD detector installed in Titan G2
%From file SSD_windowless_efficiency.xlsx
%Info is referred from FEI Appliction Note - Titan G2 with ChemiSTEM Technology A revolution in atomic analytics.
%Can also input info from other types of detector by modifying the input file. 
%Weizong Xu, April, 2015, wxu4@ncsu.edu;

eff_A = 1e-2; %initial a very small number
eff_B = 1e-2;

% if exist(filename, 'file')
%     eff_table = xlsread (filename);
    eff_table = Excel_input (filename);
if (eff_table ~=0)
    for i=1:length(eff_table)-1
        if (eff_table(i,1) > energy_A)
            %y=ax+b;
            a = (eff_table(i,2)-eff_table(i+1,2))/(eff_table(i,1)-eff_table(i+1,1));
            b = eff_table(i,2)-a*eff_table(i,1);
            eff_A = a * energy_A + b;
            break;
        end
    end
    
    for j=1:length(eff_table)-1
        if (eff_table(j,1) > energy_B)
            a = (eff_table(j,2)-eff_table(j+1,2))/(eff_table(j,1)-eff_table(j+1,1));
            b = eff_table(j,2)-a*eff_table(j,1);
            eff_B = a * energy_B + b;
            break;
        end
    end    
    

else
    uiwait(msgbox('WARNING: Unable to find input file for SSD detection efficiency, set efficiency as 1 for all X-ray energy'));
    eff_A = 1;
    eff_B = 1;
end



end

