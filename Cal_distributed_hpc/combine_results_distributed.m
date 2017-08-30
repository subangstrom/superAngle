%conbine districuted cal results
clear;close all;clc
load ('save_distributed_1.mat')
A_map_all=POS_map.A_map;
B_map_all=POS_map.B_map;
A_map_counts_all=POS_map.A_map_counts;
B_map_counts_all=POS_map.B_map_counts;
ratio_map_all=POS_map.ratio_map;
%%
Y_range=search_PosY;
X_range(1)=search_PosX(1);
for i=2:nodes_num
    load (['save_distributed_',num2str(i),'.mat'])
    if exist('x','var')
        A_map_all=[A_map_all POS_map.A_map];
        B_map_all=[B_map_all POS_map.B_map];
        A_map_counts_all=[A_map_counts_all POS_map.A_map_counts];
        B_map_counts_all=[B_map_counts_all POS_map.B_map_counts];
        ratio_map_all=[ratio_map_all POS_map.ratio_map];
    end
    
    if exist('y','var')
        A_map_all=[POS_map.A_map; A_map_all];
        B_map_all=[POS_map.B_map; B_map_all];
        A_map_counts_all=[POS_map.A_map_counts;A_map_counts_all];
        B_map_counts_all=[POS_map.B_map_counts;B_map_counts_all];
        ratio_map_all=[POS_map.ratio_map;ratio_map_all];
    end
end
if exist('x','var')
    search_PosX(1)=X_range(1);    
end

if exist('y','var')
    search_PosY(1)=Y_range(1);    
end

POS_map.A_map = A_map_all;
POS_map.B_map = B_map_all;   
POS_map.A_map_counts = A_map_counts_all;
POS_map.B_map_counts = B_map_counts_all;
POS_map.ratio_map=ratio_map_all;

chk_x=size(A_map_all,2);
chk_y=size(A_map_all,1);
if (chk_x==length(list_POSX) && chk_y==length(list_POSY))
    disp('Pass check, combined all data')
    save('POS_map_MgO_cubic_fine.mat','Detector','POS_map','chk_print','d_PosX','d_PosY','sample_para','search_PosX','search_PosY')
else
    disp('Data is not complete, please check!')
    return;
end

%% Display Position X/Y search
parameters.corr_probe=594000000/sample_para.probe_Ne; %0.1nA; check sample_para
parameters.corr_time=100/sample_para.acquire_time; %100s
parameters.corr_geo=1/input_basic.geo_real_ratio;
parameters.search_PosX=search_PosX;
parameters.d_PosX=d_PosX;
parameters.search_PosY=search_PosY;
parameters.d_PosY=d_PosY;
parameters.chk_print=2;
parameters.rot_correction=-90; %correction for titan, clockwise '-'

parameters.disp_indv_ratioMap=1; %0: not display, others: display
parameters.auto_scale_ratioMap=0; %0: manual scale, others: auto scale
parameters.Z_scale_ratioMap=[0.5,0.7];
parameters.background_level_ratioMap=0.3; %0-1 relative to the scale

parameters.disp_indv_countsMap=0; %0: not display, others: display
parameters.auto_scale_countsMap=0; %0: manual scale, others: auto scale
parameters.Z_scale_indvA_countsMap=[40000,51000];
parameters.Z_scale_indvB_countsMap=[60000,70000];
parameters.Z_scale_sumA_countsMap=[150000,180000];
parameters.Z_scale_sumB_countsMap=[250000,270000];
parameters.background_level_countsMap=0;


sample_PosXY_display2_counts_Titan( POS_map, Detector, sample_para, parameters)
