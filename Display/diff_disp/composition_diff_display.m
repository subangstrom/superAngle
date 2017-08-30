function [ ] = composition_diff_display( chk_print, comp_ratio_out, sample_para, chkX_Y )
%Display result from composition calculation
%Weizong Xu, April, 2015

if (abs(chkX_Y-2)<0.0001)
    xaxis='Tilt Y (deg)';
else
    xaxis='Tilt X (deg)';
end
[~,sym_A] = get_element_name(sample_para.EleA_num,sample_para.EleA_shell);
[~,sym_B] = get_element_name(sample_para.EleB_num,sample_para.EleB_shell);
c_select=[0.7,0.7,0;1,0,0;0,0,1;0,0.9,0;1,0,1;0,0.5,0.5;0.5,0.5,0.5;0.3,0.3,0.3;0.3,0,0.3;0.3,0.3,0;0,0.3,0.3]; %color matrix for display
Atomic_weight_A = get_element_weight(sample_para.EleA_num);
Atomic_weight_B = get_element_weight(sample_para.EleB_num);
ideal_atomic_ratio = sample_para.Atomic_ratio;
ideal_weight_ratio = sample_para.Atomic_ratio*Atomic_weight_A/Atomic_weight_B;
ideal_at_per = ideal_atomic_ratio/(1+ideal_atomic_ratio);
ideal_wt_per = ideal_weight_ratio/(1+ideal_weight_ratio);

temp_size = size (comp_ratio_out);
tot_Det_num = temp_size (2) - 2;

%convert to weight percent
comp_wt = comp_ratio_out;
comp_wt (:,2:temp_size(2)) = comp_wt (:,2:temp_size(2)) ./ (1+comp_wt (:,2:temp_size(2)));
diff_wt = comp_ratio_out;
diff_wt (:,2:temp_size(2)) =(comp_wt (:,2:temp_size(2))-ideal_wt_per)/ideal_wt_per;

temp_size2 = size(diff_wt);
y_wt_err = min(0.6,max(0.1,max(max(abs(diff_wt(:,2:temp_size2(2)))))*1.2));
y_wt_comp = max(0.1,max(max(abs(comp_wt(:,2:temp_size2(2))-ideal_wt_per)))*1.2/ideal_wt_per);


figure;
subplot(2,2,1)
hold on;
for i=1:tot_Det_num
    plot(comp_wt(:,1),comp_wt(:,i+1),'-o','color',c_select(i,:),'LineWidth',2,'MarkerSize',6);
    title(['Calculated ', sym_A, ' weight percentage x in ', sym_A, 'x', sym_B, '1-x'])
end
plot(comp_wt(:,1),comp_wt(:,tot_Det_num+2),'-o','color',[1,0,1],'LineWidth',2,'MarkerSize',6);
axis([-30 30 ideal_wt_per*(1-y_wt_comp) ideal_wt_per*(1+y_wt_comp)]);
xlabel(xaxis)
ylabel('Weight Percentage')        
grid on;
hold off;

subplot(2,2,2)
hold on;
for i=1:tot_Det_num
    plot(diff_wt(:,1),diff_wt(:,i+1),'-o','color',c_select(i,:),'LineWidth',2,'MarkerSize',6);
    title(['Calculated ', sym_A, ' weight percentage x error in ', sym_A, 'x', sym_B, '1-x'])
end
plot(diff_wt(:,1),diff_wt(:,tot_Det_num+2),'-o','color',[1,0,1],'LineWidth',2,'MarkerSize',6);
axis([-30 30 -y_wt_err y_wt_err]);
xlabel(xaxis)
ylabel('Error (a.u.)')        
grid on;
hold off;



%convert to atomic percent

comp_at = comp_ratio_out;
comp_at (:,2:temp_size(2)) = comp_at (:,2:temp_size(2)) * Atomic_weight_B/Atomic_weight_A;
comp_at (:,2:temp_size(2)) = comp_at (:,2:temp_size(2)) ./ (1+comp_at (:,2:temp_size(2)));
diff_at = comp_ratio_out;
diff_at (:,2:temp_size(2)) =(comp_at (:,2:temp_size(2))-ideal_at_per)/ideal_at_per;

temp_size2 = size(diff_at);
y_at_err = min(0.5,max(0.1,max(max(abs(diff_at(:,2:temp_size2(2)))))*1.2));
y_at_comp = max(0.1,max(max(abs(comp_at(:,2:temp_size2(2))-ideal_at_per)))*1.2/ideal_at_per);

%figure;
subplot(2,2,3)
hold on;
for i=1:tot_Det_num
    plot(comp_wt(:,1),comp_at(:,i+1),'-o','color',c_select(i,:),'LineWidth',2,'MarkerSize',6);
    title(['Calculated ', sym_A, ' atomic percentage x in ', sym_A, 'x', sym_B, '1-x'])
end
plot(comp_wt(:,1),comp_at(:,tot_Det_num+2),'-o','color',[1,0,1],'LineWidth',2,'MarkerSize',6);
axis([-30 30 ideal_at_per*(1-y_at_comp) ideal_at_per*(1+y_at_comp)]);
xlabel(xaxis)
ylabel('Atomic Percentage')        
grid on;
hold off;

subplot(2,2,4)
hold on;
for i=1:tot_Det_num
    plot(diff_at(:,1),diff_at(:,i+1),'-o','color',c_select(i,:),'LineWidth',2,'MarkerSize',6);
    title(['Calculated ', sym_A, ' atomic percentage x error in ', sym_A, 'x', sym_B, '1-x'])
end
plot(diff_at(:,1),diff_at(:,tot_Det_num+2),'-o','color',[1,0,1],'LineWidth',2,'MarkerSize',6);
axis([-30 30 -y_at_err y_at_err]);
xlabel(xaxis)
ylabel('Error (a.u.)')        
grid on;
hold off;

%comp_wt
%comp_at
diff_wt



%********************************
%set(gca,'DefaultTextFontSize',38)
c1_select=[241,196,90;190,29,44;0,144,189;137,198,55;126,47,142;162,20,47;77,190,238;0,0,0]/255;
%c1_select=[237,177,32;217,83,25;0,144,189;119,172,48;126,47,142;162,20,47;77,190,238;0,0,0]/255;
%c2_select=[245,211,131;237,139,97;81,211,255;169,214,109;198,130,213;235,90,119;139,213,243;0,0,0]/255;
%c2_select=[245,211,131;190,29,44;81,211,255;169,214,109;198,130,213;235,90,119;139,213,243;0,0,0]/255;
%c22_select=[0.7,0.7,0.35;1,0.5,0.5;0.5,0.5,1;0.45,0.9,0.45;1,0.5,1;0.25,0.5,0.5;0.25,0.25,0.25;0.15,0.15,0.15;0.3,0.15,0.3;0.3,0.3,0.15;0.15,0.3,0.3];
    %;155,30,45
c_select = [c1_select;c_select];
%c2_select = [c2_select;c22_select];
%c2_select = c_select;
%kk=1;
if (chk_print>0)
figure;
set(gca,'FontSize',18)
hold on;
plot(-30:1:30,ones(61)*0.1,'color',[0.75,0.75,0.75],'LineWidth',2);
plot(-30:1:30,-ones(61)*0.1,'color',[0.75,0.75,0.75],'LineWidth',2);
plot(-30:1:30,-ones(61)*0.0,'--','color',[0.75,0.75,0.75],'LineWidth',2);
%sym_shape=['s';'o';'^';'v';'d';'+';'*';'<';'>'];
sym_shape=['o';'o';'o';'o';'o';'o';'o';'o';'o'];


%for i=1:tot_Det_num
%    plot(diff_wt(:,1),diff_wt(:,i+1),'-o','color',c_select(i,:),'LineWidth',2,'MarkerSize',6);
%    title(['Calculated ', sym_A, ' weight percentage x error in ', sym_A, 'x', sym_B, '1-x'])
%end
%plot(diff_wt(:,1),diff_wt(:,tot_Det_num+2),'-o','color',[1,0,1],'LineWidth',2,'MarkerSize',6);
%axis([-30 30 -y_wt_err y_wt_err]);
%xlabel(xaxis)
%ylabel('Error (a.u.)')        
%grid on;
%hold off;


for i=1:tot_Det_num
    plot(diff_wt(:,1),diff_wt(:,i+1),'color',c_select(i,:),'LineWidth',3);
    %scatter(diff_wt(:,1),diff_wt(:,i+1),120,'filled',sym_shape(i,1),...
    %          'MarkerEdgeColor',c_select(i,:)*0.5,...
    %          'MarkerFaceColor',c_select(i,:), ...
    %          'LineWidth',1.5);

    title(['Difference between exp and simu in ', sym_A, ' composition (wt%)'])
end 

plot(diff_wt(:,1),diff_wt(:,tot_Det_num+2),'color',c_select(5,:),'LineWidth',5);
%scatter(diff_wt(:,1),diff_wt(:,tot_Det_num+2),120,'filled',sym_shape(i+1,1),...
%              'MarkerEdgeColor',c_select(i+1,:)*0.5,...
%              'MarkerFaceColor',c_select(i+1,:), ...
%              'LineWidth',1.5);

axis([-30 30 -y_wt_err y_wt_err]);
set(gca,'YTick',-floor(y_wt_err/0.1)*0.5*10:0.1:floor(y_wt_err/0.1)*0.5*10)
xlabel(xaxis)
ylabel('Deviation (a.u.)')        
%grid on;
box on;
%hline1 = refline([0 -0.1]);
%set(hline1,'Color',[0.75,0.75,0.75],'LineWidth',1.5)
%hline2 = refline([0 0.1]);
%set(hline2,'Color',[0.75,0.75,0.75],'LineWidth',1.5)
hold off;
if (chk_print==1)
print('wobb_Fig_out_wt_error','-depsc','-tiff')
print('wobb_Fig_out_wt_error','-dpng','-r300')
end

end





end

