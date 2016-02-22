function [diff_wt_A, diff_wt_B ] = composition_display_absolute( chk_print, comp_A_out,comp_B_out, sample_para, chkX_Y )
%Display result from composition calculation
%Weizong Xu, April, 2015

if (abs(chkX_Y-2)<0.0001)
    xaxis='Tilt Y (deg)';
else
    xaxis='Tilt X (deg)';
end
[temp, sym_A] = get_element_name(sample_para(13),sample_para(14));
[temp, sym_B] = get_element_name(sample_para(15),sample_para(16));
c_select=[0.7,0.7,0;1,0,0;0,0,1;0,0.9,0;0.5,0,0.5;0,0.5,0.5;0.5,0.5,0.5]; %color matrix for display
Atomic_weight_A = get_element_weight(sample_para(13));
Atomic_weight_B = get_element_weight(sample_para(15));
ideal_atomic_ratio = sample_para(17);
ideal_weight_ratio = sample_para(17)*Atomic_weight_A/Atomic_weight_B;
ideal_at_per_A = ideal_atomic_ratio/(1+ideal_atomic_ratio);
ideal_at_per_B = 1 - ideal_at_per_A;
ideal_wt_per_A = ideal_weight_ratio/(1+ideal_weight_ratio);
ideal_wt_per_B = 1 - ideal_wt_per_A;
temp_size = size (comp_A_out);

tot_Det_num = temp_size (2) - 2;
if (tot_Det_num>=1)

%convert to weight percent
comp_wt_A = comp_A_out;
diff_wt_A = comp_A_out;
diff_wt_A (:,2:temp_size(2)) =(comp_wt_A (:,2:temp_size(2))-ideal_wt_per_A)/ideal_wt_per_A;

temp_size2 = size(diff_wt_A);
y_wt_err_A = min(0.75,max(0.1,max(max(abs(diff_wt_A(:,2:temp_size2(2)))))*1.2));
y_wt_comp_A = max(0.1,max(max(abs(comp_wt_A(:,2:temp_size2(2))-ideal_wt_per_A)))*1.2/ideal_wt_per_A);



%convert to atomic percent

comp_at_A = comp_A_out;
comp_at_A (:,2:temp_size(2)) = comp_at_A (:,2:temp_size(2)) * Atomic_weight_B/Atomic_weight_A;
diff_at_A = comp_A_out;
diff_at_A (:,2:temp_size(2)) =(comp_at_A (:,2:temp_size(2))-ideal_at_per_A)/ideal_at_per_A;

temp_size2 = size(diff_at_A);
y_at_err_A = min(0.5,max(0.1,max(max(abs(diff_at_A(:,2:temp_size2(2)))))*1.2));
y_at_comp_A = max(0.1,max(max(abs(comp_at_A(:,2:temp_size2(2))-ideal_at_per_A)))*1.2/ideal_at_per_A);



%********************************
c1_select=[241,196,90;190,29,44;0,144,189;137,198,55;126,47,142;162,20,47;77,190,238;0,0,0]/255;
c_select = [c1_select;c_select];

if (chk_print>0)
figure;
set(gca,'FontSize',18)
hold on;
plot(-30:1:30,ones(61)*0.1,'--','color',[0.55,0.55,0.55],'LineWidth',2);
plot(-30:1:30,-ones(61)*0.1,'--','color',[0.55,0.55,0.55],'LineWidth',2);
plot(-30:1:30,-ones(61)*0.0,'--','color',[0.55,0.55,0.55],'LineWidth',2);
sym_shape=['o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o';'o'];
for i=1:tot_Det_num
    plot(diff_wt_A(:,1),diff_wt_A(:,i+1),'color',c_select(i,:),'LineWidth',2);
    scatter(diff_wt_A(:,1),diff_wt_A(:,i+1),120,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);

    title([sym_A, 'Compositional error (wt.%) of ', sym_A, ' compared with ideal Ni3Al'])
end
plot(diff_wt_A(:,1),diff_wt_A(:,tot_Det_num+2),'color',c_select(i+1,:),'LineWidth',3);
scatter(diff_wt_A(:,1),diff_wt_A(:,tot_Det_num+2),120,'filled',sym_shape(i+1,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
y_wt_err_A=0.63;
axis([-30 30 -y_wt_err_A y_wt_err_A]);
set(gca,'YTick',-floor(y_wt_err_A/0.1)*0.5:0.1:floor(y_wt_err_A/0.1)*0.5)
xlabel(xaxis)
ylabel('Deviation')        
box on;
hold off;
if (chk_print==1)
print('Fig_out_wt_error_absolute_A','-dpng','-r300')
end

end



%************for B-element display*******************

%convert to weight percent
comp_wt_B = comp_B_out;
diff_wt_B = comp_B_out;
diff_wt_B (:,2:temp_size(2)) =(comp_wt_B (:,2:temp_size(2))-ideal_wt_per_B)/ideal_wt_per_B;

temp_size2 = size(diff_wt_B);
y_wt_err_B = min(0.5,max(0.1,max(max(abs(diff_wt_B(:,2:temp_size2(2)))))*1.2));
y_wt_comp_B = max(0.1,max(max(abs(comp_wt_B(:,2:temp_size2(2))-ideal_wt_per_B)))*1.2/ideal_wt_per_B);



%convert to atomic percent

comp_at_B = comp_B_out;
comp_at_B (:,2:temp_size(2)) = comp_at_B (:,2:temp_size(2)) * Atomic_weight_B/Atomic_weight_B;
diff_at_B = comp_B_out;
diff_at_B (:,2:temp_size(2)) =(comp_at_B (:,2:temp_size(2))-ideal_at_per_B)/ideal_at_per_B;

temp_size2 = size(diff_at_B);
y_Bt_err_B = min(0.5,max(0.1,max(max(abs(diff_at_B(:,2:temp_size2(2)))))*1.2));
y_Bt_comp_B = max(0.1,max(max(abs(comp_at_B(:,2:temp_size2(2))-ideal_at_per_B)))*1.2/ideal_at_per_B);


%***************to figure*****************
if (chk_print>0)
figure;
hold on;
plot(-30:1:30,ones(61)*0.1,'--','color',[0.55,0.55,0.55],'LineWidth',2);
plot(-30:1:30,-ones(61)*0.1,'--','color',[0.55,0.55,0.55],'LineWidth',2);
plot(-30:1:30,-ones(61)*0.0,'--','color',[0.55,0.55,0.55],'LineWidth',2);
for i=1:tot_Det_num
    plot(diff_wt_B(:,1),diff_wt_B(:,i+1),'color',c_select(i,:),'LineWidth',2);
    scatter(diff_wt_B(:,1),diff_wt_B(:,i+1),120,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);

    title([sym_B, 'Compositional error (wt.%) of ', sym_B, ' compared with ideal Ni3Al'])
end
plot(diff_wt_B(:,1),diff_wt_B(:,tot_Det_num+2),'color',c_select(i+1,:),'LineWidth',3);
scatter(diff_wt_B(:,1),diff_wt_B(:,tot_Det_num+2),120,'filled',sym_shape(i+1,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
y_wt_err_B=0.63;
axis([-30 30 -y_wt_err_B y_wt_err_B]);
set(gca,'YTick',-floor(y_wt_err_B/0.1)*0.5:0.1:floor(y_wt_err_B/0.1)*0.5)
xlabel(xaxis)
ylabel('Deviation')        
box on;
hold off;
if (chk_print==1)
print('Fig_out_wt_error_Absolute_B','-dpng','-r300')
end

end









%************************end**************************

diff_wt_A

end

end

