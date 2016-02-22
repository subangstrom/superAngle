function [diff_wt] = composition_display( chk_print, comp_ratio_out, sample_para, chkX_Y )
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
ideal_at_per = ideal_atomic_ratio/(1+ideal_atomic_ratio);
ideal_wt_per = ideal_weight_ratio/(1+ideal_weight_ratio);

temp_size = size (comp_ratio_out);

tot_Det_num = temp_size (2) - 2;
if (tot_Det_num>=1)

%convert to weight percent
comp_wt = comp_ratio_out;
comp_wt (:,2:temp_size(2)) = comp_wt (:,2:temp_size(2)) ./ (1+comp_wt (:,2:temp_size(2)));
diff_wt = comp_ratio_out;
diff_wt (:,2:temp_size(2)) =(comp_wt (:,2:temp_size(2))-ideal_wt_per)/ideal_wt_per;

temp_size2 = size(diff_wt);
y_wt_err = min(0.6,max(0.1,max(max(abs(diff_wt(:,2:temp_size2(2)))))*1.2));
y_wt_comp = max(0.1,max(max(abs(comp_wt(:,2:temp_size2(2))-ideal_wt_per)))*1.2/ideal_wt_per);


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
    plot(diff_wt(:,1),diff_wt(:,i+1),'color',c_select(i,:),'LineWidth',2);
    scatter(diff_wt(:,1),diff_wt(:,i+1),120,'filled',sym_shape(i,1),...
              'MarkerEdgeColor',c_select(i,:)*0.5,...
              'MarkerFaceColor',c_select(i,:), ...
              'LineWidth',1.5);

    title([sym_A, 'Deviation between calculated composition (Al wt.%) and ideal value of Ni3Al (Ratio method)'])
end
plot(diff_wt(:,1),diff_wt(:,tot_Det_num+2),'color',c_select(i+1,:),'LineWidth',3);
scatter(diff_wt(:,1),diff_wt(:,tot_Det_num+2),120,'filled',sym_shape(i+1,1),...
              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
              'MarkerFaceColor',c_select(i+1,:), ...
              'LineWidth',1.5);
y_wt_err=0.63;
axis([-30 30 -y_wt_err y_wt_err]);
set(gca,'YTick',-floor(y_wt_err/0.1)*0.5:0.1:floor(y_wt_err/0.1)*0.5)
xlabel(xaxis)
ylabel('Deviation')        
box on;
hold off;
if (chk_print==1)
print('Fig_out_wt_error','-dpng','-r300')
end

end



%comp_wt
%comp_at
diff_wt

end

end

