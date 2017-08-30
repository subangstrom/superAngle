% control_file_generate for distributed run on hpc cluster
clear;close all;clc
%***************************************************************************************************************
input_basic.model_filename='model_connect_cubic_195_205_190.mat';
input_basic.save_filename='save_MgO.mat';
input_basic.core_num=4; %number of cores
input_basic.dAngle=0.5; %deg
input_basic.detector_filename='detector.xlsx';
input_basic.specimen_filename='specimen_MgO_particle.xlsx';
input_basic.holder_filename='holder_FEI_LB.xlsx';
%input_basic.t_chk=1; % 1 --> constant thickness t, other value --> constant spot during tilt, i.e. t will change
input_basic.Z_rotation=12+90; %degree
input_basic.geo_real_ratio=1; % thickness calculation for nanoparticles, i.e. set sample size. %nm in unit
input_basic.t_chk=0; %=0 use model defined thickness;
input_basic.chk_Shadow=2; 
%>0 and <>2) consider simple holder shadowing (circular contours and fully cutoff)
%2) consider holder shadowing by FEI low background double tilt HIVis holder
%<=0 value) No shadowing effect

%***************************************************************************************************************
input_SPOT.run=0; %1 for run, others skip
input_SPOT.ether_start_p=[0,0,0];
input_SPOT.chk_2nd=-1; %1 consider 2nd interception (25-100x slow), other not consider
input_SPOT.TiltX=0;
input_SPOT.TiltY=0;

%***************************************************************************************************************
input_POS_search.run=1; %1 for run, others skip
input_POS_search.search_PosX=[-120,120];
input_POS_search.search_PosY=[-120,120];
input_POS_search.d_PosX = 1;
input_POS_search.d_PosY = 1;
input_POS_search.chk_print = 2; %1 output high quality image to file, >0 display but not output, <0 not display
input_POS_search.chk_2nd=-1; %1 consider 2nd interception (25-100x slow), other not consider

%***************************************************************************************************************
%% setup number of nodes for distribution calculation
nodes_num=8;
if (nodes_num>=2)
    list_POSX=input_POS_search.search_PosX(1):input_POS_search.d_PosX:input_POS_search.search_PosX(2);
    %segment along y direction
    list_POSY=input_POS_search.search_PosY(1):input_POS_search.d_PosY:input_POS_search.search_PosY(2);
    if (length(list_POSX)>length(list_POSY))
        dx_num=ceil(length(list_POSX)/nodes_num);
        for i=1:nodes_num
            x(1)=list_POSX(dx_num*(i-1)+1);
            x(2)=list_POSX(min(dx_num*i,length(list_POSX)));
            disp([num2str(x),'  ',num2str(x(2)-x(1)+1)])
            input_POS_search.search_PosX=[x(1),x(2)];
            input_basic.save_filename=['save_distributed_',num2str(i),'.mat'];
            save (['control_',num2str(i),'.mat'])
        end
    else
        dy_num=ceil(length(list_POSY)/nodes_num);
        for i=1:nodes_num
            y(1)=list_POSY(dy_num*(i-1)+1);
            y(2)=list_POSY(min(dy_num*i,length(list_POSY)));
            disp([num2str(y),'  ',num2str(y(2)-y(1)+1)])
            input_POS_search.search_PosY=[y(1),y(2)];
            input_basic.save_filename=['save_distributed_',num2str(i),'.mat'];
            save (['control_',num2str(i),'.mat'])
        end
    end
else
    save ('control.mat')
end

%% test if all nodes works
if (input_SPOT.run==0 && input_POS_search.run==0)
    for i=1:nodes_num
        run_filename=['main_script_V2_nanoParticle_control_hpc',num2str(i),'.m'];
        run(run_filename);
    end
end