% control_file_generate for single run
clear;close all;clc
%***************************************************************************************************************
input_basic.model_filename='model_connect_rod_100_rect.mat';
input_basic.save_filename='save_single.mat';
input_basic.core_num=4; %number of cores
input_basic.dAngle=2; %deg
input_basic.detector_filename='detector.xlsx';
input_basic.specimen_filename='specimen_wedge_ideal.xlsx';
input_basic.holder_filename='holder_FEI_LB.xlsx';
%input_basic.t_chk=1; % 1 --> constant thickness t, other value --> constant spot during tilt, i.e. t will change
input_basic.Z_rotation=0; %degree
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
input_POS_search.search_PosX=[-150,150];
input_POS_search.search_PosY=[-50,50];
input_POS_search.d_PosX = 10;
input_POS_search.d_PosY = 10;
input_POS_search.chk_print = 2; %1 output high quality image to file, >0 display but not output, <0 not display
input_POS_search.chk_2nd=-1; %1 consider 2nd interception (25-100x slow), other not consider

%***************************************************************************************************************
save ('control.mat')