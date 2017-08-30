function [ ] = line_display_Counts_paper( line_search_Result, chk_print, chkXY, Detector, exp_file, search_Deg, sample_para )
%Display 1D plot for correction coefficient due to X-ray absorption and/or holder shadowing
%The result will also be compared with the experiment data.
%Weizong Xu, wxu4@ncsu.edu, March 2015

A_line=line_search_Result.A_line;
B_line=line_search_Result.B_line;
A_Abso=line_search_Result.A_Abso;
B_Abso=line_search_Result.B_Abso;
Absrp_line=line_search_Result.Absrp_line;

A_lineall_counts(:,1)=A_line(:,1,1);
A_lineall_counts(:,2)=A_line(:,2,1);
B_lineall_counts(:,1)=B_line(:,1,1);
B_lineall_counts(:,2)=B_line(:,2,1);
tot_Det_num=Detector.tot_Det_num;
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

sym_A = get_element_name(sample_para.EleA_num,sample_para.EleA_shell);
sym_B = get_element_name(sample_para.EleB_num,sample_para.EleB_shell);



figure;
subplot(2,4,1)

hold on;
for i=1:tot_Det_num
    Plot2Dim(A_Abso(:,:,i), ['Absorption percentage ', sym_A], 'p', xaxis, 'a.u.');
    axis([-search_Deg search_Deg 0 1])
end
hold off;

subplot(2,4,2)
%figure;
hold on;
for i=1:tot_Det_num
    Plot2Dim(B_Abso(:,:,i), ['Absorption percentage ', sym_B], 'p', xaxis, 'a.u.');
    axis([-search_Deg search_Deg 0 1])
end
hold off;

subplot(2,4,3)
%figure;
hold on;
for i=1:tot_Det_num
    Plot2Dim(A_line(:,:,i), [sym_A, ' total counts'], 'p', xaxis, 'a.u.');
end
scatter(A_lineall_counts(:,1),A_lineall_counts(:,2),'filled',...
              'MarkerEdgeColor',[0 0.5 1],...
              'MarkerFaceColor',[0 0.5 1]);
hold off;

subplot(2,4,4)
%figure;
hold on;
for i=1:tot_Det_num
    Plot2Dim(B_line(:,:,i), [sym_B, ' total counts'], 'p', xaxis, 'a.u.');
end
scatter(B_lineall_counts(:,1),B_lineall_counts(:,2),'filled',...
              'MarkerEdgeColor',[0 0.5 1],...
              'MarkerFaceColor',[0 0.5 1]);
hold off;


%Data input from xlsx file e.g. 'exp_B3A_wedge.xlsx'; 
%if file exist -->  tot_Det_num_in > 0
[ A_exp, A_exp_norm, max_A_exp, B_exp, B_exp_norm, max_B_exp, ratio_exp, ratio_exp_all, tot_Det_num_in] = read_exp_data( exp_file );
if (tot_Det_num_in == 0)
    figure;
    hold on;
    for i=1:tot_Det_num
        Plot2Dim(Absrp_line(:,:,i), 'Absorption correction coefficient', 'p', xaxis, 'a.u.');
        axis([-search_Deg search_Deg 0.0 1])
    end
    scatter(Absrp_lineall(:,1),Absrp_lineall(:,2),'b','filled');
    hold off;
    
    return;
end

A_exp_all_norm(:,1)=A_exp_norm(:,1,1);
A_exp_all_norm(:,2)=A_exp_norm(:,2,1);
B_exp_all_norm(:,1)=B_exp_norm(:,1,1);
B_exp_all_norm(:,2)=B_exp_norm(:,2,1);
for i=2:tot_Det_num_in
    A_exp_all_norm(:,2)=A_exp_all_norm(:,2)+A_exp_norm(:,2,i);
    B_exp_all_norm(:,2)= B_exp_all_norm(:,2)+B_exp_norm(:,2,i);
end


    
    max_A=max(max(A_line(:,2,:)));
    max_B=max(max(B_line(:,2,:)));
    for i=1:tot_Det_num
        A_line_norm(:,1,i) = A_line(:,1,i);
        A_line_norm(:,2,i) = A_line(:,2,i)/max_A;

        B_line_norm(:,1,i) = B_line(:,1,i);
        B_line_norm(:,2,i) = B_line(:,2,i)/max_B;
        %A_line_sum_norm = A_line_sum_norm + A_line_norm(:,2,i); 
    end
    
    
    
    c_select=[0.7,0.7,0;1,0,0;0,0,1;0,0.9,0;1,0,1;0,0.5,0.5;0.5,0.5,0.5;0.3,0.3,0.3;0.3,0,0.3;0.3,0.3,0;0,0.3,0.3]; %color matrix for display
    subplot(2,4,5)
    %figure;
    hold on;
    for i=1:tot_Det_num
        Plot2Dim(A_line_norm(:,:,i), [sym_A, ' normalized counts'], 'p', xaxis, 'a.u.');
    end
    
    for i=1:tot_Det_num_in
        scatter(A_exp_norm(:,1,i),A_exp_norm(:,2,i),'filled',...
              'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',c_select(i,:));
    end
    axis([-search_Deg search_Deg 0 1.1])
    hold off;

    subplot(2,4,6)
    %figure;
    hold on;
    for i=1:tot_Det_num
       Plot2Dim(B_line_norm(:,:,i), [sym_B, ' normalized counts'], 'p', xaxis, 'a.u.');
    end
    
    for i=1:tot_Det_num_in
               scatter(B_exp_norm(:,1,i),B_exp_norm(:,2,i),'filled',...
              'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',c_select(i,:));

    end    
    
    axis([-search_Deg search_Deg 0 1.1])
    hold off;

    subplot(2,4,7)
    %figure;
    hold on;
    for i=1:tot_Det_num
        Plot2Dim(A_line_norm(:,:,i), [sym_A, ' normalized counts'], 'p', xaxis, 'a.u.');
    end
    scatter(A_lineall_counts(:,1),A_lineall_counts(:,2)/max_A,'filled',...
              'MarkerEdgeColor',[0 0.5 1],...
              'MarkerFaceColor',[0 0.5 1]);
          
    for i=1:tot_Det_num_in
        scatter(A_exp_norm(:,1,i),A_exp_norm(:,2,i),'filled',...
              'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',c_select(i,:));
    end
    
        scatter(A_exp_all_norm(:,1),A_exp_all_norm(:,2),'filled',...
              'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',c_select(i+1,:));
    
    axis([-search_Deg search_Deg 0 4.2])
    hold off;

    subplot(2,4,8)
    %figure;
    hold on;
    for i=1:tot_Det_num
       Plot2Dim(B_line_norm(:,:,i), [sym_B, ' normalized counts'], 'p', xaxis, 'a.u.');
    end
    scatter(B_lineall_counts(:,1),B_lineall_counts(:,2)/max_B,'filled',...
              'MarkerEdgeColor',[0 0.5 1],...
              'MarkerFaceColor',[0 0.5 1]);
          
    for i=1:tot_Det_num_in
               scatter(B_exp_norm(:,1,i),B_exp_norm(:,2,i),'filled',...
              'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',c_select(i,:));

    end
    
        scatter(B_exp_all_norm(:,1),B_exp_all_norm(:,2),'filled',...
              'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',c_select(i+1,:));
    
    axis([-search_Deg search_Deg 0 4.2])
    hold off;
    
    
    
    figure;
    subplot(2,2,1)
    hold on;
    for i=1:tot_Det_num
        Plot2Dim(Absrp_line(:,:,i), 'Correction coefficient due to absorption/shadowing', 'p', xaxis, 'a.u.');
    end
    scatter(Absrp_lineall(:,1),Absrp_lineall(:,2),'filled',...
              'MarkerEdgeColor',[0 0 1],...
              'MarkerFaceColor',[0 0.5 1]);
    hold off;

    %factor ~=1/3.7;%dependes on ion crossection, detector efficiency, probe current...
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
    
%    figure;
    subplot(2,2,2)
    hold on;
    for i=1:tot_Det_num
        Plot2Dim(Absrp_line(:,:,i), ['Calculated counts ratio of ' sym_A, '/', sym_B], 'p', xaxis, 'a.u.');
    end
    scatter(Absrp_lineall(:,1),Absrp_lineall(:,2),'filled',...
              'MarkerEdgeColor',[0 0 1],...
              'MarkerFaceColor',[0 0.5 1]);
    hold off;
    
    range_exp_max = max(max(max(ratio_exp(:,2,:))),max(ratio_exp_all(:,2)))*1.1;
    subplot(2,2,3)
    %figure;
    hold on;
    for i=1:tot_Det_num_in
        plot(ratio_exp(:,1,i),ratio_exp(:,2,i),'-o','color',c_select(i,:),'LineWidth',2,'MarkerSize',8);
        title(['Experimenal counts ratio of ', sym_A, '/', sym_B])
    end
    plot(ratio_exp_all(:,1),ratio_exp_all(:,2),'-o','color',[1,0,1],'LineWidth',2,'MarkerSize',8);
    axis([-search_Deg search_Deg 0 range_exp_max]);
    xlabel(xaxis)
    ylabel('a.u.')        
    grid on;
    hold off;

    range_simu_max = max(max(max(Absrp_line(:,2,:))),max(Absrp_lineall(:,2)))*1.1;
    range_max = max(range_exp_max,range_simu_max);
    subplot(2,2,4)
    %figure;
    hold on;
    for i=1:tot_Det_num
       Plot2Dim(Absrp_line(:,:,i), 'Comparison between Exp and Cal counts ratio', 'p', xaxis, 'a.u.');
    end
    
    for i=1:tot_Det_num_in
       plot(ratio_exp(:,1,i),ratio_exp(:,2,i),'-o','color',c_select(i,:),'LineWidth',2,'MarkerSize',8);
    end
    scatter(Absrp_lineall(:,1),Absrp_lineall(:,2),'filled',...
              'MarkerEdgeColor',[0 0 1],...
              'MarkerFaceColor',[0 0.5 1]);
    plot(ratio_exp_all(:,1),ratio_exp_all(:,2),'-o','color',[1,0,1],'LineWidth',2,'MarkerSize',8);
    axis([-search_Deg search_Deg 0 range_max])
    hold off;

    
    
%********************************************************   
c1_select=[241,196,90;190,29,44;0,144,189;137,198,55;126,47,142;162,20,47;77,190,238;0,0,0]/255;
%c1_select=[237,177,32;217,83,25;0,144,189;119,172,48;126,47,142;162,20,47;77,190,238;0,0,0]/255;
%c2_select=[245,211,131;237,139,97;81,211,255;169,214,109;198,130,213;235,90,119;139,213,243;0,0,0]/255;
%c2_select=[245,211,131;190,29,44;81,211,255;169,214,109;198,130,213;235,90,119;139,213,243;0,0,0]/255;
%c22_select=[0.7,0.7,0.35;1,0.5,0.5;0.5,0.5,1;0.45,0.9,0.45;1,0.5,1;0.25,0.5,0.5;0.25,0.25,0.25;0.15,0.15,0.15;0.3,0.15,0.3;0.3,0.3,0.15;0.15,0.3,0.3];
    %;155,30,45
c_select = [c1_select;c_select];
%c2_select = [c2_select;c22_select];
c2_select = c_select;
%color_all=[171,119,61]/255; %for total counts
color_all=c_select(5,:);
%sym_shape=['s';'o';'^';'v';'d';'+';'*';'<';'>'];
sym_shape=['o';'o';'o';'o';'o';'o';'o';'o';'o'];
%kk=1;
%chk_print = -1;

if (chk_print>0)
figure;
%subplot(2,4,1)
set(gca,'FontSize',18)
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(A_line(:,:,i), ['Effective solid angle of ', sym_A, ' X-ray'], 3, xaxis, 'a.u.',c_select(i,:),1);
    axis([-search_Deg search_Deg 0 0.2])
end
box on;
hold off;
if (chk_print==1)
print('Fig_out_Abs_perc_Al','-depsc','-tiff');
print('Fig_out_Abs_perc_Al','-dpng','-r300')
end


%subplot(2,4,2)
figure;
set(gca,'FontSize',18)
hold on;
for i=1:tot_Det_num
    Plot2Dim_factor(B_line(:,:,i), ['Effective solid angle of ', sym_B, 'X-ray'], 3, xaxis, 'a.u.',c_select(i,:),1);
    axis([-search_Deg search_Deg 0 0.2])
end
box on;
hold off;
if (chk_print==1)
print('Fig_out_Abs_perc_Ni','-depsc','-tiff')
print('Fig_out_Abs_perc_Ni','-dpng','-r300')
end

%subplot(2,4,3)
%figure;
%hold on;
%for i=1:tot_Det_num
%    Plot2Dim(A_line(:,:,i), [sym_A, ' total counts'], 'p', xaxis, 'a.u.');
%end
%scatter(A_lineall_counts(:,1),A_lineall_counts(:,2),'filled',...
%              'MarkerEdgeColor',[0 0.5 1],...
%              'MarkerFaceColor',[0 0.5 1]);
%hold off;

%subplot(2,4,4)
%figure;
%hold on;
%for i=1:tot_Det_num
%    Plot2Dim(B_line(:,:,i), [sym_B, ' total counts'], 'p', xaxis, 'a.u.');
%end
%scatter(B_lineall_counts(:,1),B_lineall_counts(:,2),'filled',...
%              'MarkerEdgeColor',[0 0.5 1],...
%              'MarkerFaceColor',[0 0.5 1]);
%hold off;

%    subplot(2,4,5)

%probe_current = 0.14;%nA
%convert_factor_A = probe_current*acquire_time*sample_para(5)*90.33;%Al %only valid for constant thickness condition now, need construction later
%convert_factor_B = probe_current*acquire_time*sample_para(5)*346.28;%Ni

[ convert_factor_A,convert_factor_B, ~,~ ] = absolute_scale_factor( sample_para );

A_line_count(:,1,:) = A_line(:,1,:);
A_line_count(:,2,:) = A_line(:,2,:)*convert_factor_A;
B_line_count(:,1,:) = B_line(:,1,:);
B_line_count(:,2,:) = B_line(:,2,:)*convert_factor_B;
A_line_count_all(:,1) = A_lineall_counts(:,1,:);
A_line_count_all(:,2) = A_lineall_counts(:,2,:)*convert_factor_A;
B_line_count_all(:,1) = B_lineall_counts(:,1,:);
B_line_count_all(:,2) = B_lineall_counts(:,2,:)*convert_factor_B;

    figure;
    set(gca,'FontSize',18)
    hold on;
    
    for i=1:tot_Det_num
        Plot2Dim_factor(A_line_count(:,:,i), ['Compare with experiment data - ', sym_A], 3, xaxis, 'Counts a.u.',c_select(i,:),1);%0.97
    end
    
    for i=1:tot_Det_num_in
        add_error_counts(c_select(i,:), A_exp(:,:,i),1);   
        scatter(A_exp(:,1,i),A_exp(:,2,i),120,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);
    end
    %axis([-search_Deg search_Deg 0 1.1])
    hold off;
    box on;
    
if (chk_print==1)
print('Fig_out_counts_compare_Al','-depsc','-tiff');
print('Fig_out_counts_compare_Al','-dpng','-r300')
end


    %subplot(2,4,6)
    figure;
    set(gca,'FontSize',18)
    hold on;
    for i=1:tot_Det_num
       Plot2Dim_factor(B_line_count(:,:,i), ['Compare with experimental data - ', sym_B], 3, xaxis, 'Counts a.u.', c_select(i,:),1);%1
    end
    
    for i=1:tot_Det_num_in
            add_error_counts(c_select(i,:), B_exp(:,:,i),1);
            scatter(B_exp(:,1,i),B_exp(:,2,i),120, 'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);
    end    
    
    %axis([-search_Deg search_Deg 0 1.1])
    hold off;
    box on;
if (chk_print==1)
print('Fig_out_counts_compare_Ni','-depsc','-tiff');
print('Fig_out_counts_compare_Ni','-dpng','-r300')
end


    exp_all_A(:,1)=A_exp(:,1,1);
    exp_all_A(:,2)=A_exp(:,2,1);
    exp_all_B(:,1)=B_exp(:,1,1);
    exp_all_B(:,2)=B_exp(:,2,1);
    for i=2:tot_Det_num
        exp_all_A(:,2)=exp_all_A(:,2)+A_exp(:,2,i);
        exp_all_B(:,2)=exp_all_B(:,2)+B_exp(:,2,i);
    end


    figure;
    set(gca,'FontSize',18)
    hold on;
    for i=1:tot_Det_num
       Plot2Dim_factor(A_line_count(:,:,i), ['Compare with experimental data - ', sym_A], 3, xaxis, 'Counts a.u.', c_select(i,:),1);%1
    end
Plot2Dim_factor(A_line_count_all, ['Compare with experimental data - ', sym_A], 5, xaxis, 'Counts a.u.', c_select(i+1,:),1);%1

    for i=1:tot_Det_num_in
            add_error_counts(c_select(i,:), A_exp(:,:,i),1); 
            scatter(A_exp(:,1,i),A_exp(:,2,i),120, 'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);
    end
    add_error_counts(color_all, A_exp,2);
    scatter(exp_all_A(:,1),exp_all_A(:,2),120, 'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
    
    %axis([-search_Deg search_Deg 0 1.1])
    hold off;
    box on;
if (chk_print==1)
print('Fig_out_counts_compare_all_Al','-depsc','-tiff');
print('Fig_out_counts_compare_all_Al','-dpng','-r300')
end


    figure;
    set(gca,'FontSize',18)
    hold on;
    Plot2Dim_factor(A_line_count_all, ['Compare with total data - ', sym_A], 5, xaxis, 'Counts a.u.', c_select(i+1,:),1);%1
    add_error_counts(color_all, A_exp,2);
    scatter(exp_all_A(:,1),exp_all_A(:,2),120, 'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
    
    %axis([-search_Deg search_Deg 0 1.1])
    hold off;
    box on;
if (chk_print==1)
print('Fig_out_counts_compare_total_Al','-depsc','-tiff');
print('Fig_out_counts_compare_total_Al','-dpng','-r300')
end


    figure;
    set(gca,'FontSize',18)
    hold on;
    for i=1:tot_Det_num
       add_error_counts(c_select(i,:), B_exp(:,:,i),1); 
       Plot2Dim_factor(B_line_count(:,:,i), ['Compare with experimental data - ', sym_B], 3, xaxis, 'Counts a.u.', c_select(i,:),1);%1
    end
Plot2Dim_factor(B_line_count_all, ['Compare with experimental data - ', sym_B], 5, xaxis, 'Counts a.u.', c_select(i+1,:),1);%1

    for i=1:tot_Det_num_in
            add_error_counts(color_all, B_exp,2);
            scatter(B_exp(:,1,i),B_exp(:,2,i),120, 'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);
    end    
            add_error_counts(color_all, B_exp,2);
            scatter(exp_all_B(:,1),exp_all_B(:,2),120, 'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
    
    %axis([-search_Deg search_Deg 0 1.1])
    hold off;
    box on;
if (chk_print==1)
print('Fig_out_counts_compare_all_Ni','-depsc','-tiff');
print('Fig_out_counts_compare_all_Ni','-dpng','-r300')
end


    figure;
    set(gca,'FontSize',18)
    set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);
    hold on;
    Plot2Dim_factor(A_line_count_all, ['Compare with total data - ', sym_A], 5, xaxis, 'Counts a.u.', color_all,1);%1
    add_error_counts(color_all, A_exp,2);
    scatter(exp_all_A(:,1),exp_all_A(:,2),120, 'filled',sym_shape(i,1),...
              'MarkerEdgeColor',color_all*0.5,...
              'MarkerFaceColor',color_all, ...
              'LineWidth',1.5);
    y_max_Al=max(exp_all_A(:,2))*1.1;
    y_min_Al=min(A_line_count_all(:,2))*0.9;
    axis([-search_Deg search_Deg y_min_Al y_max_Al])
    hold off;
    box on;
if (chk_print==1)
set(gcf, 'PaperPositionMode', 'auto') 
print('Fig_out_counts_compare_total_Al','-depsc','-tiff');
print('Fig_out_counts_compare_total_Al','-dpng','-r300')
end

    figure;
    set(gca,'FontSize',18)
    %set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.32]);
    set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);
    hold on;
    %color_all=[171,119,61]/255;
    Plot2Dim_factor(B_line_count_all, ['Compare with total data - ', sym_B], 5, xaxis, 'Counts a.u.', color_all,1);%1
    add_error_counts(color_all, B_exp,2);
    scatter(exp_all_B(:,1),exp_all_B(:,2),120, 'filled',sym_shape(i,1),...
              'MarkerEdgeColor',color_all*0.5,...
              'MarkerFaceColor',color_all, ...
              'LineWidth',1.5);
    y_max_Ni=max(exp_all_B(:,2))*1.1;
    y_min_Ni=min(B_line_count_all(:,2))*0.9;
    axis([-search_Deg search_Deg y_min_Ni y_max_Ni])
    hold off;
    box on;
if (chk_print==1)
set(gcf, 'PaperPositionMode', 'auto') 
print('Fig_out_counts_compare_total_Ni','-depsc','-tiff');
print('Fig_out_counts_compare_total_Ni','-dpng','-r300')
end

%     subplot(2,4,7)
%     %figure;
%     hold on;
%     for i=1:tot_Det_num
%         Plot2Dim(A_line_norm(:,:,i), [sym_A, ' normalized counts'], 'p', xaxis, 'a.u.');
%     end
%     scatter(A_lineall_counts(:,1),A_lineall_counts(:,2)/max_A,'filled',...
%               'MarkerEdgeColor',[0 0.5 1],...
%               'MarkerFaceColor',[0 0.5 1]);
%           
%     for i=1:tot_Det_num_in
%         scatter(A_exp_norm(:,1,i),A_exp_norm(:,2,i),'filled',...
%               'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',c_select(i,:));
%     end
%     
%         scatter(A_exp_all_norm(:,1),A_exp_all_norm(:,2),'filled',...
%               'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',c_select(i+1,:));
%     
%     axis([-search_Deg search_Deg 0 4.2])
%     hold off;
% 
%     subplot(2,4,8)
%     %figure;
%     hold on;
%     for i=1:tot_Det_num
%        Plot2Dim(B_line_norm(:,:,i), [sym_B, ' normalized counts'], 'p', xaxis, 'a.u.');
%     end
%     scatter(B_lineall_counts(:,1),B_lineall_counts(:,2)/max_B,'filled',...
%               'MarkerEdgeColor',[0 0.5 1],...
%               'MarkerFaceColor',[0 0.5 1]);
%           
%     for i=1:tot_Det_num_in
%                scatter(B_exp_norm(:,1,i),B_exp_norm(:,2,i),'filled',...
%               'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',c_select(i,:));
% 
%     end
%     
%         scatter(B_exp_all_norm(:,1),B_exp_all_norm(:,2),'filled',...
%               'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',c_select(i+1,:));
%     
%     axis([-search_Deg search_Deg 0 4.2])
%     hold off;
    
    
    
%     figure;
%     subplot(2,2,1)
%     hold on;
%     for i=1:tot_Det_num
%         Plot2Dim(Absrp_line(:,:,i), 'Correction coefficient due to absorption/shadowing', 'p', xaxis, 'a.u.');
%     end
%     scatter(Absrp_lineall(:,1),Absrp_lineall(:,2),'filled',...
%               'MarkerEdgeColor',[0 0 1],...
%               'MarkerFaceColor',[0 0.5 1]);
%     hold off;

    
    figure;
    set(gca,'FontSize',18)
    %subplot(2,2,2)
    hold on;
    for i=1:tot_Det_num
        Plot2Dim_factor(Absrp_line(:,:,i), ['Calculated counts ratio of ' sym_A, '/', sym_B], 3, xaxis, 'a.u.',c_select(i,:),1);
    end
            Plot2Dim_factor(Absrp_lineall, ['Calculated counts ratio of ' sym_A, '/', sym_B], 5, xaxis, 'a.u.',c2_select(i+1,:),1);
legend('detector 1','detector 2','detector 3','detector 4','all detectors','Location','northoutside','Orientation','horizontal')

%    scatter(Absrp_lineall(:,1),Absrp_lineall(:,2),'filled',...
%              'MarkerEdgeColor',[0 0 1],...
%              'MarkerFaceColor',[0 0.5 1]);
    hold off;
    box on;
if (chk_print==1)
print('Fig_out_ratio_cal','-depsc','-tiff');
print('Fig_out_ratio_1cal','-dpng','-r300')
end

%    subplot(2,2,3)
    figure;
    set(gca,'FontSize',18)
    hold on;
    for i=1:tot_Det_num_in
        plot(ratio_exp(:,1,i),ratio_exp(:,2,i),'color',c_select(i,:),'LineWidth',2);
        add_error_counts(c_select(i,:), A_exp(:,:,i),3,B_exp(:,:,i));
        scatter(ratio_exp(:,1,i),ratio_exp(:,2,i),120,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);
        title(['Counts ratio of ', sym_A, '/', sym_B, ' from experiment'])
    end
    plot(ratio_exp_all(:,1),ratio_exp_all(:,2),'color',c_select(i+1,:),'LineWidth',3);
    add_error_counts(c_select(i+1,:), A_exp,4,B_exp);
    scatter(ratio_exp_all(:,1),ratio_exp_all(:,2),120,'filled',sym_shape(i+1,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
    axis([-search_Deg search_Deg 0 range_exp_max]);
    xlabel(xaxis)
    ylabel('a.u.')        
    %grid on;
    box on;
    hold off;
if (chk_print==1)
print('Fig_out_ratio_exp','-depsc','-tiff');
print('Fig_out_ratio_2exp','-dpng','-r300')
end

    %subplot(2,2,4)
    figure;
    set(gca,'FontSize',18)
    hold on;
    for i=1:tot_Det_num
       Plot2Dim_factor(Absrp_line(:,:,i), 'Comparison between Exp and Cal counts ratio', 3, xaxis, 'a.u.',c2_select(i,:),1);
    end

    for i=1:tot_Det_num_in
       %plot(ratio_exp(:,1,i),ratio_exp(:,2,i),'-o','color',c_select(i,:),'LineWidth',2,'MarkerSize',8);
%       plot(ratio_exp(:,1,i),ratio_exp(:,2,i),'color',c_select(i,:),'LineWidth',1);
        add_error_counts(c_select(i,:), A_exp(:,:,i),3,B_exp(:,:,i));       
        scatter(ratio_exp(:,1,i),ratio_exp(:,2,i),120,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);
    end
    Plot2Dim_factor(Absrp_lineall, ['Calculated counts ratio of ' sym_A, '/', sym_B], 5, xaxis, 'a.u.',c2_select(i+1,:),1);

    %scatter(Absrp_lineall(:,1),Absrp_lineall(:,2),'filled',...
    %          'MarkerEdgeColor',[0 0 1],...
    %          'MarkerFaceColor',[0 0.5 1]);
    %plot(ratio_exp_all(:,1),ratio_exp_all(:,2),'-o','color',[1,0,1],'LineWidth',2,'MarkerSize',8);
%    plot(ratio_exp_all(:,1),ratio_exp_all(:,2),'color',c_select(i+1,:),'LineWidth',2);
    add_error_counts(c_select(i+1,:), A_exp,4,B_exp);
    scatter(ratio_exp_all(:,1),ratio_exp_all(:,2),120,'filled',sym_shape(i+1,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
    axis([-search_Deg search_Deg 0 range_max])
    box on;
    hold off;
%legend('detctor1','detector2','detector3','detector4','all detectors','2','3','4','5','6','Location','northoutside','Orientation','horizontal')
    
if (chk_print==1)
    print('Fig_out_ratio_compare','-depsc','-tiff');
    print('Fig_out_ratio_3compare','-dpng','-r300');
end
end

end
