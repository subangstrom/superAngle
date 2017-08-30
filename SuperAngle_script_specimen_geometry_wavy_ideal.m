%% Calculate Abosrption and Holder Shadowing Correction Factors for Super-X EDS detectors
% This script is particulary for non-ideal or distorted specimen
% Version 2.0beta (script version)
% Considering
% (a) Specimen absorption
% (b) Shadowing from holder frame (Currently FEI low background holder)
% (c) Selective filtering from holder carrier (Be stage)
% (d) Shadowing from Softlock securing clip
% (e) Shadowing from Supporting grid (circle or 1mmx2mm slot)
% (f) 3D based model for more complicated specimen geometry
% To import model, please run code accordingly in 3D_model_build folder
% Consider real holder geometry
% New data structure and processing engine for Matlab 2014b later version
% Support for multiple detectors, default for Super-X EDS detectors.
% Support for parallel computing
% GUI version is under construction
% Modeled and Coded by Weizong Xu
% Email to wxu4@ncsu.edu for any information and suggestions. 
% August 2017

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
load ('3D_model_sinwave.mat');
%% rotate model if needed,
axis=[0,0,1];theta=180;center=[0,100,0];
[ p ] = rotate_omni_center( axis, theta, center, model_sample.p');
model_sample.p=p';
[ p ] = rotate_omni_center( axis, theta, center, model_ether.p');
model_ether.p=p';
%% **********Setup Detector***************
dAngle = 0.5; %Accuracy for angular intergration, 0.5 is fairly good for speed, 0.2 is better for accuracy
[ Detector ] = Detector_input( 'detector.xlsx', dAngle );
% %*****Spurious Xray calculation ***********
%Be cautious!!! Set to zero if not sure
[ SpuriousX ] = SpuriousXray_setup( Detector.tot_Det_num, 0, 0 );
%% **********Setup specimen************
t_chk = 1; % 1 --> constant thickness t input from sample parameter, 0- use model defined thickness(const spot) (nm), 10- match thickness with model defined thickness (search spot), other value --> constant spot during tilt, i.e. t will change
[ sample_para ] = Specimen_setup( 'specimen_wedge_ideal.xlsx', t_chk );
%[ sample_para ] = Specimen_setup( 'specimen_wedge_real.xlsx', t_chk );
sample_para.Thickness=50;
sample_para.Z_rotation=0; %degree
[ model_sample, model_ether ]= tilt_model_precal( 0, 0, sample_para.Z_rotation, model_sample, model_ether);
model_sample.geo_real_ratio=1/25; % thickness calculation for nanoparticles, i.e. set sample size. %nm in unit
sample_para.chk_2nd=-1; %1 consider 2nd interception (25-100x slow), other not consider
%**********Setup holder**************
chk_Shadow = 2;
%>0 and <>2) consider simple holder shadowing (circular contours and fully cutoff)
%2) consider holder shadowing by FEI low background double tilt HIVis holder
%<=0 value) No shadowing effect
[ holder_para ] = Holder_setup( chk_Shadow, 'holder_FEI_LB.xlsx', sample_para );
%[ holder_para ] = Holder_setup( chk_Shadow, 'holder_FEI_LB_demo_sample.xlsx', sample_para );

%% +++++Invididual Spot Check (@X Y tilt)
sample_para_tmp=sample_para;
sample_para_tmp.Deviation_angle_X=0;
sample_para_tmp.Deviation_angle_Y=0;
sample_para_tmp.chk_2nd=-1; %1 consider 2nd interception (25-100x slow), other not consider

ether_start_p=[50,100,0]; %starting point of the model relative to ether (if sample shape varies, this could be useful)
tic
[ output_Point1 ] = single_spot_3D_display(0, 0, Detector, sample_para_tmp, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
toc

%% +++++Do TiltX TiltY search from -search_Deg to +search_Deg
sample_para.Deviation_angle_X=0;
sample_para.Deviation_angle_Y=0;
search_Deg=30;
d_Deg = 1;
exp_file = 'exp_Ni3Al_wedge.xlsx';
chkXY = 1; %2-Search Y, others search X
chk_print = 2; %1 output high quality image to file, others>0 display but not output, <0 not display

sample_para.chk_2nd=-1; %1 consider 2nd interception (25-100x slow), other not consider
ether_start_p=[-50,50,0]; %starting point of the model relative to ether (if sample shape varies, this could be useful)
chkXY = 1;
[ line_search_Result_X ] = line_search_3D_fast( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
chkXY = 2;
[ line_search_Result_Y ] = line_search_3D_fast( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);

%ideal case for comparision
[ line_search_Result_std_X ] = line_search( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX);
chkXY = 2;
[ line_search_Result_std_Y ] = line_search( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX);

%% display difference
chk_print = 1; %1 output high quality image to file, others>0 display but not output, <0 not display
chkXY=1; %1-X 2-Y
line_search_Result_std=line_search_Result_std_X; % assume geometry is real, this is the standard
line_search_Result=line_search_Result_X; % assume flat, deviated from real
diff_display( chk_print, line_search_Result, line_search_Result_std, chkXY, Detector, search_Deg, sample_para);
%% display some data
line_display_Counts( line_search_Result_std_X, chk_print, chkXY, Detector, exp_file, search_Deg, sample_para );
