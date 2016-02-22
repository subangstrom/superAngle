function [ comp_A_out,comp_B_out] = composition_cal_absolute( exp_file, tot_Det_num1, sample_para, holder_para, holder_frame_para, angle_search, SpuriousX, chkX_Y )
%Calculate composition of AxBy based on k_ideal and correction coefficient
%Weizong Xu, April, 2015


if exist(exp_file, 'file')
[ A_exp, A_exp_norm, max_A_exp, B_exp, B_exp_norm, max_B_exp, ratio_exp, ratio_exp_all, tot_Det_num] = read_exp_data( exp_file );

exp_all_A(:,1)=A_exp(:,1,1);
exp_all_A(:,2)=A_exp(:,2,1);
exp_all_B(:,1)=B_exp(:,1,1);
exp_all_B(:,2)=B_exp(:,2,1);
for i=2:tot_Det_num
   exp_all_A(:,2)=exp_all_A(:,2)+A_exp(:,2,i);
   exp_all_B(:,2)=exp_all_B(:,2)+B_exp(:,2,i);
end



if (tot_Det_num ==0)
    comp_A_out=0;
    comp_B_out=0;
    return;
end

if (tot_Det_num > tot_Det_num1)
    uiwait(msgbox('Error!! Number of detectors in Exp is larger than that in Simulation!!'));
    comp_A_out=0;
    comp_B_out=0;
    return;
end

if (tot_Det_num < tot_Det_num1)
    uiwait(msgbox('WARNING!! Number of detectors in Exp is less than that in Simu!'));
end

cal_chk=sample_para(12); 

if (cal_chk ~= 1)

    comp_A_out = zeros(length(exp_all_A), tot_Det_num+2);
    comp_B_out = zeros(length(exp_all_B), tot_Det_num+2);
    comp_A_out(:,1)=exp_all_A(:,1);
    comp_B_out(:,1)=exp_all_B(:,1);
    for i=1:length(A_exp)
    
        if (chkX_Y == 2) % 2) along Y tilt
            TiltY=A_exp(i,1);
            TiltX=sample_para(3);
        else
            TiltX=A_exp(i,1); % other value, along X tilt
            TiltY=sample_para(4);
        end
        
        [tempSP, point_out] = single_spot(TiltX, TiltY, tot_Det_num, sample_para, holder_para, holder_frame_para, angle_search, SpuriousX);
        omega_A = point_out(:,2);
        omega_B = point_out(:,3);
        [ tempA, tempB, convert_factor_A, convert_factor_B ] = absolute_scale_factor( sample_para );
       
        
        for j=1:tot_Det_num
            comp_A_out (i,j+1)=A_exp(i,2,j)/convert_factor_A/omega_A(j);
            comp_B_out (i,j+1)=B_exp(i,2,j)/convert_factor_B/omega_B(j);
        end
            comp_A_out (i,tot_Det_num+2)=exp_all_A(i,2)/convert_factor_A/omega_A(tot_Det_num+1);
            comp_B_out (i,tot_Det_num+2)=exp_all_B(i,2)/convert_factor_B/omega_B(tot_Det_num+1);

    end
    
else
%
    
end


else
    comp_A_out=0;
    comp_B_out=0;

end

end





