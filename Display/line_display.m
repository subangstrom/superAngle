function [ ] = line_display( line_search_Result, chk_print, chkXY, Detector, exp_file, search_Deg, sample_para )
%Display 1D plot for correction coefficient due to X-ray absorption and/or holder shadowing
%The result will also be compared with the experiment data.
%Weizong Xu, wxu4@ncsu.edu, March 2015

A_line=line_search_Result.A_line;
B_line=line_search_Result.B_line;
A_Abso=line_search_Result.A_Abso;
B_Abso=line_search_Result.B_Abso;
Absrp_line=line_search_Result.Absrp_line;

A_line_all(:,1)=A_line(:,1,1);
A_line_all(:,2)=A_line(:,2,1);
B_line_all(:,1)=B_line(:,1,1);
B_line_all(:,2)=B_line(:,2,1);
tot_Det_num=Detector.tot_Det_num;
for i=2:tot_Det_num
    A_line_all(:,2)=A_line_all(:,2)+A_line(:,2,i);
    B_line_all(:,2)= B_line_all(:,2)+B_line(:,2,i);
end
Absrp_lineall(:,1)=A_line_all(:,1);
Absrp_lineall(:,2)=A_line_all(:,2)./B_line_all(:,2);

%Display results

if (abs(chkXY(1)-2)<0.0001)
    xaxis='Tilt Y (deg)';
else
    xaxis='Tilt X (deg)';
end

sym_A = get_element_name(sample_para.EleA_num,sample_para.EleA_shell);
sym_B = get_element_name(sample_para.EleB_num,sample_para.EleB_shell);

k_AB_ideal = sample_para.k_AB_ideal; % ideal k ratio in for composition in wt%
atomic_ratio_AB = sample_para.Atomic_ratio;
Atomic_weight_A = get_element_weight(sample_para.EleA_num);
Atomic_weight_B = get_element_weight(sample_para.EleB_num);
k_AB_ideal_atomic = k_AB_ideal/Atomic_weight_A*Atomic_weight_B;
Int_ratio_factor= atomic_ratio_AB/k_AB_ideal_atomic;

for i=1:tot_Det_num
    Absrp_line(:,2,i) = Absrp_line(:,2,i)*Int_ratio_factor;
end
Absrp_lineall(:,2) = Absrp_lineall(:,2)*Int_ratio_factor;

    
    
%********************************************************   
c1_select=[241,196,90;190,29,44;0,144,189;137,198,55;126,47,142;162,20,47;77,190,238;0,0,0]/255;
c_select=[0.7,0.7,0;1,0,0;0,0,1;0,0.9,0;1,0,1;0,0.5,0.5;0.5,0.5,0.5;0.3,0.3,0.3;0.3,0,0.3;0.3,0.3,0;0,0.3,0.3]; %color matrix for display
c_select = [c1_select;c_select];
c2_select = c_select;

figure;
subplot(2,2,1)
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(A_Abso(:,:,i), ['Absorption percentage ', sym_A, ' X-ray'], 3, xaxis, '%',c_select(i,:),100);
    xlim([-search_Deg search_Deg]);
end
hold off;

subplot(2,2,2)
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(B_Abso(:,:,i), ['Absorption percentage ', sym_B, ' X-ray'], 3, xaxis, '%',c_select(i,:),100);
    xlim([-search_Deg search_Deg]);
end
hold off;

subplot(2,2,3)
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(A_line(:,:,i), ['Effective solid angle of ', sym_A, ' X-ray'], 3, xaxis, 'Omega (Sr)',c_select(i,:),1);
    xlim([-search_Deg search_Deg]);
end
hold on;
    Plot2Dim_factor(A_line_all, ['Effective solid angle of ', sym_A, ' X-ray'], 3, xaxis,  'Omega (Sr)', c_select(i+1,:),1);
hold off;

subplot(2,2,4)
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(B_line(:,:,i), ['Effective solid angle of ', sym_B, 'X-ray'], 3, xaxis,  'Omega (Sr)', c_select(i,:),1);
    xlim([-search_Deg search_Deg]);
end
hold on;
Plot2Dim_factor(B_line_all, ['Effective solid angle of ', sym_B, ' X-ray'], 3, xaxis,  'Omega (Sr)', c_select(i+1,:),1);
hold off;
if (chk_print==1)
    print('Line_display_solid_angle','-dpng','-r300')
end

figure;
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(Absrp_line(:,:,i), ['Calculated counts ratio of ' sym_A, '/', sym_B], 3, xaxis, [sym_A, '/', sym_B],c_select(i,:),1);
end
    Plot2Dim_factor(Absrp_lineall, ['Calculated counts ratio of ' sym_A, '/', sym_B], 5, xaxis, [sym_A, '/', sym_B],c2_select(i+1,:),1);
xlim([-search_Deg search_Deg]);
for i=1:tot_Det_num
    legend_str{i}=['detector ', num2str(i)];
end
legend_str{i+1}='Sum';
legend(legend_str,'Location','northoutside','Orientation','horizontal');
hold off;
box on;
if (chk_print==1)
    print('Line_display_ratio_no_exp','-dpng','-r300')
end

end