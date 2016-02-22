function [ ] = line_display_Counts( chk_print, chkXY, tot_Det_num, exp_file, search_Deg, A_line, B_line, A_Abso, B_Abso, Absrp_line, sample_para )
%Display 1D plot for correction coefficient due to X-ray absorption and/or holder shadowing
%The result will also be compared with the experiment data.
%Weizong Xu, wxu4@ncsu.edu, April 2015

A_lineall_counts(:,1)=A_line(:,1,1);
A_lineall_counts(:,2)=A_line(:,2,1);
B_lineall_counts(:,1)=B_line(:,1,1);
B_lineall_counts(:,2)=B_line(:,2,1);
for i=2:tot_Det_num
    A_lineall_counts(:,2)=A_lineall_counts(:,2)+A_line(:,2,i);
    B_lineall_counts(:,2)= B_lineall_counts(:,2)+B_line(:,2,i);
end
Absrp_lineall(:,1)=A_lineall_counts(:,1);
Absrp_lineall(:,2)=A_lineall_counts(:,2)./B_lineall_counts(:,2);

%Display results

if (abs(chkXY(1)-2)<0.0001)
    xaxis='Tilt Y (deg)';
else
    xaxis='Tilt X (deg)';
end

sym_A = get_element_name(sample_para(13),sample_para(14));
sym_B = get_element_name(sample_para(15),sample_para(16));

k_AB_ideal = sample_para(20); % ideal k ratio in for composition in wt%
atomic_ratio_AB = sample_para(17);
Atomic_weight_A = get_element_weight(sample_para(13));
Atomic_weight_B = get_element_weight(sample_para(15));
k_AB_ideal_atomic = k_AB_ideal/Atomic_weight_A*Atomic_weight_B;
Int_ratio_factor= atomic_ratio_AB/k_AB_ideal_atomic;
for i=1:tot_Det_num
    Absrp_line(:,2,i) = Absrp_line(:,2,i)*Int_ratio_factor;
end
    Absrp_lineall(:,2) = Absrp_lineall(:,2)*Int_ratio_factor;

    
[ convert_factor_A,convert_factor_B, tempA, tempB ] = absolute_scale_factor(sample_para);
A_line_count(:,1,:) = A_line(:,1,:);
A_line_count(:,2,:) = A_line(:,2,:)*convert_factor_A;
B_line_count(:,1,:) = B_line(:,1,:);
B_line_count(:,2,:) = B_line(:,2,:)*convert_factor_B;
A_line_count_all(:,1) = A_lineall_counts(:,1,:);
A_line_count_all(:,2) = A_lineall_counts(:,2,:)*convert_factor_A;
B_line_count_all(:,1) = B_lineall_counts(:,1,:);
B_line_count_all(:,2) = B_lineall_counts(:,2,:)*convert_factor_B;    
    
%********************************************************   
c1_select=[241,196,90;190,29,44;0,144,189;137,198,55;126,47,142;162,20,47;77,190,238;0,0,0]/255;
c_select=[0.7,0.7,0;1,0,0;0,0,1;0,0.9,0;1,0,1;0,0.5,0.5;0.5,0.5,0.5;0.3,0.3,0.3;0.3,0,0.3;0.3,0.3,0;0,0.3,0.3]; %color matrix for display
c_select = [c1_select;c_select];
c2_select = c_select;
color_all=c_select(5,:);
sym_shape=['o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o'];

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

if (sample_para(22)>0 && sample_para(23)>0) %must have value of tau and e current
subplot(2,2,3)
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(A_line_count(:,:,i), ['Calculated counts of ', sym_A, ' X-ray'], 3, xaxis, 'Counts',c_select(i,:),1);
    xlim([-search_Deg search_Deg]);
end
hold on;
    Plot2Dim_factor(A_line_count_all, ['Calculated counts of ', sym_A, ' X-ray'], 3, xaxis,  'Counts', c_select(i+1,:),1);
hold off;

subplot(2,2,4)
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(B_line_count(:,:,i), ['Calculated counts of ', sym_B, 'X-ray'], 3, xaxis, 'Counts',c_select(i,:),1);
    xlim([-search_Deg search_Deg]);
end
hold on;
Plot2Dim_factor(B_line_count_all, ['Calculated counts of ', sym_B, ' X-ray'], 3, xaxis,  'Counts', c_select(i+1,:),1);
hold off;
end %end of tau and current check

figure;
%set(gca,'FontSize',18)
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
    
%Data input from xlsx file e.g. 'exp_B3A_wedge.xlsx'; 
%if file exist -->  tot_Det_num_in > 0
[ A_exp, A_exp_norm, max_A_exp, B_exp, B_exp_norm, max_B_exp, ratio_exp, ratio_exp_all, tot_Det_num_in] = read_exp_data( exp_file );
if (tot_Det_num_in == 0 || sample_para(22)<=0 || sample_para(23)<=0)
    return; % no further experimental comparison, return
end


%************Compare with exp and simu*****************
    
figure;
    subplot(2,2,1)
    hold on;
    
    for i=1:tot_Det_num
        Plot2Dim_factor(A_line_count(:,:,i), ['Compare with experiment data - ', sym_A], 3, xaxis, 'Counts a.u.',c_select(i,:),1);%0.97
    end
    
    for i=1:tot_Det_num_in
        add_error_counts(c_select(i,:), A_exp(:,:,i),1);        
        scatter(A_exp(:,1,i),A_exp(:,2,i),80,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);

    end
    xlim([-search_Deg search_Deg]);
    hold off;
    box on;



    subplot(2,2,2)
    hold on;
    for i=1:tot_Det_num
       Plot2Dim_factor(B_line_count(:,:,i), ['Compare with experimental data - ', sym_B], 3, xaxis, 'Counts a.u.', c_select(i,:),1);%1
    end
    
    for i=1:tot_Det_num_in
        add_error_counts(c_select(i,:), B_exp(:,:,i),1);        
        scatter(B_exp(:,1,i),B_exp(:,2,i),80,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);

    end    
    
    xlim([-search_Deg search_Deg]);
    hold off;
    box on;


    exp_all_A(:,1)=A_exp(:,1,1);
    exp_all_A(:,2)=A_exp(:,2,1);
    exp_all_B(:,1)=B_exp(:,1,1);
    exp_all_B(:,2)=B_exp(:,2,1);
    for i=2:tot_Det_num
        exp_all_A(:,2)=exp_all_A(:,2)+A_exp(:,2,i);
        exp_all_B(:,2)=exp_all_B(:,2)+B_exp(:,2,i);
    end

    subplot(2,2,3)
    hold on;
    Plot2Dim_factor(A_line_count_all, ['Compare with total data - ', sym_A], 5, xaxis, 'Counts a.u.', color_all,1);%1
    add_error_counts(color_all, A_exp,2);
    scatter(exp_all_A(:,1),exp_all_A(:,2),80,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',color_all*0.5,...
              'MarkerFaceColor',color_all, ...
              'LineWidth',1.5);

    y_max_Al=max(exp_all_A(:,2))*1.1;
    y_min_Al=min(A_line_count_all(:,2))*0.9;
    axis([-search_Deg search_Deg y_min_Al y_max_Al])
    hold off;
    box on;
    
    subplot(2,2,4)
    hold on;
    Plot2Dim_factor(B_line_count_all, ['Compare with total data - ', sym_B], 5, xaxis, 'Counts a.u.', color_all,1);%1
    add_error_counts(color_all, B_exp,2);
    scatter(exp_all_B(:,1),exp_all_B(:,2),80,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',color_all*0.5,...
              'MarkerFaceColor',color_all, ...
              'LineWidth',1.5);
    
    y_max_Ni=max(exp_all_B(:,2))*1.1;
    y_min_Ni=min(B_line_count_all(:,2))*0.9;
    axis([-search_Deg search_Deg y_min_Ni y_max_Ni])
    hold off;
    box on;



    figure;
    hold on;
    for i=1:tot_Det_num
       Plot2Dim_factor(Absrp_line(:,:,i), 'Comparison between Exp and Cal counts ratio', 3, xaxis, 'a.u.',c2_select(i,:),1);
    end
    Plot2Dim_factor(Absrp_lineall, ['Calculated counts ratio of ' sym_A, '/', sym_B], 5, xaxis, 'a.u.',c2_select(i+1,:),1);
    legend(legend_str,'Location','northoutside','Orientation','horizontal');
    for i=1:tot_Det_num_in
       add_error_counts(c_select(i,:), A_exp(:,:,i),3,B_exp(:,:,i));
       scatter(ratio_exp(:,1,i),ratio_exp(:,2,i),80,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);

    end

    add_error_counts(c_select(i+1,:), A_exp,4,B_exp);
    scatter(ratio_exp_all(:,1),ratio_exp_all(:,2),80,'filled',sym_shape(i+1,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
    range_exp_max = max(max(max(ratio_exp(:,2,:))),max(ratio_exp_all(:,2)))*1.1;
    range_simu_max = max(max(max(Absrp_line(:,2,:))),max(Absrp_lineall(:,2)))*1.1;
    range_max = max(range_exp_max,range_simu_max);
    axis([-search_Deg search_Deg 0 range_max])

    box on;
    hold off;


end
