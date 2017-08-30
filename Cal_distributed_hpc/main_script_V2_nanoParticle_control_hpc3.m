%% Calculate Abosrption and Holder Shadowing Correction Factors for Super-X EDS detectors

%=======Coordinate========
%   Titan G2 ChemiSTEM
%           Y+
%            |
%       D1   |     D2  
%     135deg |   45deg
%           (Z)-----X+
%            
%       D4         D3
%     225deg     315deg
%            
%=========================

clear;close all;clc;
load ('control_3.mat') %will load input_basic
disp('Node #3 is running...')
load (input_basic.model_filename);
%% prepare for parallel running
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool(input_basic.core_num);
end
disp(['Note: CPU cores for this run is ',num2str(poolobj.NumWorkers),'.'])
%% **********Setup Detector***************
dAngle = input_basic.dAngle; %Accuracy for angular intergration, 0.5 is fairly good for speed, 0.2 is better for accuracy
[ Detector ] = Detector_input( input_basic.detector_filename, dAngle );
% %*****Spurious Xray calculation ***********
%Be cautious!!! Set to zero if not sure
[ SpuriousX ] = SpuriousXray_setup( Detector.tot_Det_num, 0, 0 );
%% **********Setup specimen************
t_chk = input_basic.t_chk; % 1 --> constant thickness t, other value --> constant spot during tilt, i.e. t will change
[ sample_para ] = Specimen_setup( input_basic.specimen_filename, t_chk );
%[ sample_para ] = Specimen_setup( 'specimen_wedge_real.xlsx', t_chk );
sample_para.Z_rotation=input_basic.Z_rotation; %degree
[ model_sample, model_ether ]= tilt_model_precal( 0, 0, sample_para.Z_rotation, model_sample, model_ether);
model_sample.geo_real_ratio=input_basic.geo_real_ratio; % thickness calculation for nanoparticles, i.e. set sample size. %nm in unit
sample_para.t_chk=input_basic.t_chk; %use model defined thickness;
%**********Setup holder**************
chk_Shadow = input_basic.chk_Shadow;
%>0 and <>2) consider simple holder shadowing (circular contours and fully cutoff)
%2) consider holder shadowing by FEI low background double tilt HIVis holder
%<=0 value) No shadowing effect
[ holder_para ] = Holder_setup( chk_Shadow, input_basic.holder_filename, sample_para );
%[ holder_para ] = Holder_setup( chk_Shadow, 'holder_FEI_LB_demo_sample.xlsx', sample_para );

%% +++++Invididual Spot Check (@X Y tilt)
if (input_SPOT.run==1)
sample_para.chk_2nd=input_SPOT.chk_2nd; %1 consider 2nd interception (25-100x slow), other not consider
ether_start_p=input_SPOT.ether_start_p; %starting point of the model relative to ether (if sample shape varies, this could be useful)
TiltX=input_SPOT.TiltX;
TiltY=input_SPOT.TiltY;
tic
[ output_Point11 ] = single_spot_3D_display(TiltX, TiltY, Detector, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
toc
% tic
% [ output_Point1 ] = single_spot_3D(0, 0, Detector, sample_para_tmp, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
% toc
% tic
% [ output_Point2 ] = single_spot(0, 0, Detector, sample_para_tmp, holder_para, SpuriousX); % counts*1e3;
% toc
end
%% +++++Do Position X/Y search
if (input_POS_search.run==1)
search_PosX=input_POS_search.search_PosX;
search_PosY=input_POS_search.search_PosY;
d_PosX = input_POS_search.d_PosX;
d_PosY = input_POS_search.d_PosY;
chk_print = input_POS_search.chk_print; %1 output high quality image to file, others>0 display but not output, <0 not display
sample_para.chk_2nd=input_POS_search.chk_2nd; %1 consider 2nd interception (25-100x slow), other not consider
tic
[ POS_map ] = sample_POS_search_3D_fast( Detector, search_PosX, search_PosY, d_PosX, d_PosY, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect);
toc
sample_PosXY_display2_counts( POS_map, chk_print, Detector, sample_para, search_PosX, d_PosX, search_PosY, d_PosY )
%save ('POS_map_rod1.mat', 'POS_map', 'chk_print', 'Detector', 'sample_para', 'search_PosX', 'd_PosX', 'search_PosY', 'd_PosY');
end
%% +++++Do TiltX TiltY search from -search_Deg to +search_Deg
% sample_para.Deviation_angle_X=0;
% sample_para.Deviation_angle_Y=0;
% search_Deg=30;
% d_Deg = 1;
% exp_file = 'exp_Ni3Al_wedge.xlsx';
% chkXY = 1; %2-Search Y, others search X
% chk_print = 2; %1 output high quality image to file, others>0 display but not output, <0 not display
% sample_para.chk_2nd=-1; %1 consider 2nd interception (25-100x slow), other not consider
% tic
% ether_start_p=[0,0,0]; %starting point of the model relative to ether (if sample shape varies, this could be useful)
% [ line_search_Result1 ] = line_search_3D_fast( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
% toc %275.7(0.5dA30deg1deg) 1347.45(0.2dA30deg1deg) 868(older)
% line_display_Counts( line_search_Result1, chk_print, chkXY, Detector, exp_file, search_Deg, sample_para );
% %save ('line_search_result.mat','line_search_Result1', 'chk_print', 'chkXY', 'Detector', 'exp_file', 'search_Deg', 'sample_para', 'ether_start_p')
% tic
% [ line_search_Result2 ] = line_search( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX);
% toc %6.56, 5.35(0.5dA30deg1deg) 34.42(0.2dA30deg1deg)
% line_display_Counts( line_search_Result2, chk_print, chkXY, Detector, exp_file, search_Deg, sample_para );
% 
% %speed ratio: 275.7/6.56=42.0   1347.45/34.42=39.14
% %% +++++Do 2D tiltX-Y full range search [-search_Deg , +search_Deg], must know composition
% search_Deg_2D=30;
% d_Deg_2D = 5;
% chk_print = 2;
% 
% tic
% ether_start_p=[0.15,0.15,0.15];
% [ Tilt_map1 ] = XY_search_3D( Detector, search_Deg_2D, d_Deg_2D, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
% XY_display2_solidAngle( Tilt_map1, chk_print, Detector, sample_para, search_Deg_2D, d_Deg_2D );
% XY_display2_counts( Tilt_map1, chk_print, Detector, sample_para, search_Deg_2D, d_Deg_2D );
% toc
% tic
% [ Tilt_map2 ] = XY_search( Detector, search_Deg_2D, d_Deg_2D, sample_para, holder_para, SpuriousX);
% XY_display2_solidAngle( Tilt_map2, chk_print, Detector, sample_para, search_Deg_2D, d_Deg_2D );
% XY_display2_counts( Tilt_map2, chk_print, Detector, sample_para, search_Deg_2D, d_Deg_2D );
% toc
%% save file and shut down parpool
clear model_connect model_ether model_sample %reduce output file size
save (input_basic.save_filename);
delete(poolobj);