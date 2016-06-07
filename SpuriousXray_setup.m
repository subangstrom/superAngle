function [ SpuriousX ] = SpuriousXray_setup( tot_Det_num, ele_A, ele_B )
%Adding effect from Spurious Xray
%Not very scientific in current model, need more study
%Be cautious!!! Set zero if not sure

for i=1:tot_Det_num
    
    SpuriousX(1,i)=0;%ele_A; %add arbitrary counts e.g. 500
    SpuriousX(2,i)=0;%ele_B; %add arbitrary counts e.g. 500

end

end

