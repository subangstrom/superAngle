%% load data from Display_combined_result_TiltXY_v2_range_for_figure.m
clear;
load ('out_E14.mat');
load ('out_E22.mat');
load ('out_E18.mat');
load ('out_E18_grid.mat');

%%
c_select=[241,196,90;190,29,44;0,144,189;137,198,55;126,47,142;162,20,47;77,190,238;0,0,0]/255;
c_select2=[190,29,44;0,144,189;209,152,18;73,106,30;162,20,47;77,190,238;241,196,90;126,47,142;0,0,0]/255;
color_all=c_select(5,:);
    
%% figure 1 for holder_grid_shadowing
figure
%det 1

ha20=plot(x_E18{1},y_E18{2},'color',[0,0,0],'LineWidth',1.5);
hold on;
ha21=plot(x_E18_grid{1},y_E18_grid{3},'--','color',[0,0,0],'LineWidth',1.5);
hold on;

ha0=plot(x_E18{1},y_E18{1},'color',c_select(1,:),'LineWidth',1.5);
hold on;
ha1=plot(x_E18_grid{1},y_E18_grid{1},'--','color',c_select(1,:),'LineWidth',1.5);
hold on;
hs1=fill([x_E18{1};flipud(x_E18{1})],[y_E18_grid{1};flipud(y_E18{1})],c_select(1,:));
set(hs1,'EdgeColor','None');
alpha(hs1,.1)

%det 3
hb0=plot(x_E18{1},y_E18{2},'color',c_select(3,:),'LineWidth',2);
hold on;
hb1=plot(x_E18_grid{1},y_E18_grid{2},'--','color',c_select(3,:),'LineWidth',1.5);
hold on;
hs2=fill([x_E18{1};flipud(x_E18{1})],[y_E18_grid{2};flipud(y_E18{2})],c_select(3,:));
set(hs2,'EdgeColor','None');
alpha(hs2,.05)


%det 1+3
hc0=plot(x_E18{1},y_E18{3},'color',c_select(4,:),'LineWidth',2);
hold on;
hc1=plot(x_E18_grid{1},y_E18_grid{3},'--','color',c_select(4,:),'LineWidth',1.5);
hold on;
hs3=fill([x_E18{1};flipud(x_E18{1})],[y_E18_grid{3};flipud(y_E18{3})],c_select(4,:));
set(hs3,'EdgeColor','None');
alpha(hs3,.08)

%det sum
hd0=plot(x_E18{1},y_E18{4},'color',c_select(5,:),'LineWidth',4);
hold on;
hd1=plot(x_E18_grid{1},y_E18_grid{4},'--','color',c_select(5,:),'LineWidth',2.5);
hold on;
hs4=fill([x_E18{1};flipud(x_E18{1})],[y_E18_grid{4};flipud(y_E18{4})],c_select(5,:));
set(hs4,'EdgeColor','None');
alpha(hs4,.2)

xlim([min(x_E14{1}), max(x_E14{1})])
ylim([0,50])
legend([ha0 hb0 hc0 hd0 ha20 ha21],{'Det 1','Det 3','Det 1+3','Sum','No grid shadow','Grid shadow'},'Location','northwest')
xlabel('Tilt angle (deg)'),
ylabel('Error (%)')

% print('Search_error_XY_grid_map','-depsc','-tiff');
% print('Search_error_XY_grid_map','-dpng','-r300');

%% figure 2 for elevation angle
figure
%det 1
a21=plot(x_E18{1},y_E18{1},':','color',[0,0,0],'LineWidth',1.5);
hold on;
a22=plot(x_E18{1},y_E18{2},'color',[0,0,0],'LineWidth',1.5);
hold on;
a23=plot(x_E18{1},y_E18{3},'--','color',[0,0,0],'LineWidth',1.5);
hold on;
a1=plot(x_E14{1},y_E14{1},':','color',c_select(1,:),'LineWidth',1.5);
hold on;
a2=plot(x_E18{1},y_E18{1},'color',c_select(1,:),'LineWidth',1.5);
hold on;
a3=plot(x_E22{1},y_E22{1},'--','color',c_select(1,:),'LineWidth',1.5);
hold on;
s1=fill([x_E14{1};flipud(x_E14{1})],[y_E14{1};flipud(y_E22{1})],c_select(1,:));
set(s1,'EdgeColor','None');
alpha(s1,.1)

%det 3
b1=plot(x_E14{1},y_E14{2},':','color',c_select(3,:),'LineWidth',1.5);
hold on;
b2=plot(x_E18{1},y_E18{2},'color',c_select(3,:),'LineWidth',2);
hold on;
b3=plot(x_E22{1},y_E22{2},'--','color',c_select(3,:),'LineWidth',1.5);
hold on;
s2=fill([x_E14{1};flipud(x_E14{1})],[y_E14{2};flipud(y_E22{2})],c_select(3,:));
set(s2,'EdgeColor','None');
alpha(s2,.05)


%det 1+3
c1=plot(x_E14{1},y_E14{3},':','color',c_select(4,:),'LineWidth',1.5);
hold on;
c2=plot(x_E18{1},y_E18{3},'color',c_select(4,:),'LineWidth',2);
hold on;
c3=plot(x_E22{1},y_E22{3},'--','color',c_select(4,:),'LineWidth',1.5);
hold on;
s3=fill([x_E14{1};flipud(x_E14{1})],[y_E14{3};flipud(y_E22{3})],c_select(4,:));
set(s3,'EdgeColor','None');
alpha(s3,.08)

%det sum
d1=plot(x_E14{1},y_E14{4},':','color',c_select(5,:),'LineWidth',2.5);
hold on;
d2=plot(x_E18{1},y_E18{4},'color',c_select(5,:),'LineWidth',4);
hold on;
d3=plot(x_E22{1},y_E22{4},'--','color',c_select(5,:),'LineWidth',2.5);
hold on;
s4=fill([x_E14{1};flipud(x_E14{1})],[y_E14{4};flipud(y_E22{4})],c_select(5,:));
set(s4,'EdgeColor','None');
alpha(s4,.2)

xlim([min(x_E14{1}), max(x_E14{1})])
ylim([0,50])
legend([a2 b2 c2 d2 a21 a22 a23],{'Det 1','Det 3','Det 1+3','Sum','\theta_E=14','\theta_E=18','\theta_E=22'},'Location','northwest')
xlabel('Tilt angle (deg)'),
ylabel('Error (%)')

% print('Search_error_XY_thetaE_map','-depsc','-tiff');
% print('Search_error_XY_thetaE_map','-dpng','-r300');
