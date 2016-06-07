%% Calculate Abosrption and Holder Shadowing Correction Factors for Super-X EDS detectors
% Version 1.30s (script version)
% Considering
% (a) Specimen absorption
% (b) Shadowing from holder frame (Currently FEI low background holder)
% (c) Selective filtering from holder carrier (Be stage)
% (d) Shadowing from Softlock securing clip
% (e) Shadowing from Supporting grid (circle or 1mmx2mm slot)
% New data structure and processing engine for Matlab 2014b later version
% Support for multiple detectors, default for Super-X EDS detectors.
% Support for parallel computing.
% For accurate predictions, please check the sample geometry properly.
% Modeled and Coded by Weizong Xu
% Email to wxu4@ncsu.edu for any information and suggestions. 
% June 2016

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

clear;
close all;
clc;

%**********Setup Detector***************
dAngle = 0.5; %Accuracy for angular intergration, 0.5 is fairly good for speed, 0.2 is better for accuracy
[ Detector ] = Detector_input( 'detector.xlsx', dAngle );
% %*****Spurious Xray calculation ***********
%Be cautious!!! Set to zero if not sure
ele_A =0; %for Al
ele_B =0; %for Ni
[ SpuriousX ] = SpuriousXray_setup( Detector.tot_Det_num, ele_A, ele_B );
%% **********Setup specimen************
t_chk = 1; % 1 --> constant thickness t, other value --> constant spot during tilt, i.e. t will change
[ sample_para ] = Specimen_setup( 'specimen_Ni3Al_demo.xlsx', t_chk );
%**********Setup holder**************
chk_Shadow = 2;
%>0 and <>2) consider simple holder shadowing (circular contours and fully cutoff)
%2) consider holder shadowing by FEI low background double tilt HIVis holder
%<=0 value) No shadowing effect
[ holder_para ] = Holder_setup( chk_Shadow, 'holder_FEI_LB_demo_sample.xlsx', sample_para );

%% +++++Invididual Spot Check (@X Y tilt)
[ output_Point ] = single_spot(0, 0, Detector, sample_para, holder_para, SpuriousX); % counts*1e3;
%[ output_Point ] = single_spot_counts(0, 0, tot_Det_num, sample_para, holder_para, holder_frame_para, angle_search, SpuriousX); % counts*1e3;
%calculate composition
exp.A_counts=12691;exp.B_counts=83422;
[comp_ratio_weight_Spot, comp_ratio_atomic_spot] = composition_cal_single_spot( exp, Detector, sample_para, holder_para, SpuriousX);
%% +++++Do TiltX TiltY search from -search_Deg to +search_Deg
search_Deg=30;
d_Deg = 1;
exp_file = 'Exp_data_Ni3Al_demo_X.xlsx';
chkXY = 1; %2-Search Y, others search X
chk_print = 2; %1 output high quality image to file, others>0 display but not output, <0 not display
[ line_search_Result ] = line_search( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX);
line_display_Counts( line_search_Result, chk_print, chkXY, Detector, exp_file, search_Deg, sample_para );
%% +++++Composition calculation from experimental input
if (sample_para.cal_chk==0)
[comp_ratio_out] = composition_cal( exp_file, Detector, sample_para, holder_para, SpuriousX, chkXY );
[diff_wt]=composition_display( chk_print, comp_ratio_out, sample_para, chkXY );
[comp_A_out, comp_B_out] = composition_cal_absolute( exp_file, Detector, sample_para, holder_para, SpuriousX, chkXY );
[diff_wt_A, diff_wt_B]=composition_display_absolute( chk_print, comp_A_out, comp_B_out, sample_para, chkXY );
end
%% +++++Do 2D tiltX-Y full range search [-search_Deg , +search_Deg], must know composition
search_Deg_2D=30;
d_Deg_2D = 2;
chk_print = 2;
[ Tilt_map ] = XY_search( Detector, search_Deg_2D, d_Deg_2D, sample_para, holder_para, SpuriousX);
XY_display2_solidAngle( Tilt_map, chk_print, Detector, sample_para, search_Deg_2D, d_Deg_2D );
XY_display2_counts( Tilt_map, chk_print, Detector, sample_para, search_Deg_2D, d_Deg_2D );
%% +++++Do 2D position shift X-Y full range search [-search_Range , +search_Range], must know composition
search_Range_2D=0.6; %in mm
d_Range_2D = 0.1;
chk_print = 2;
[ Shift_map ] = positionXY_search( Detector, search_Range_2D, d_Range_2D, sample_para, holder_para, SpuriousX);
PositionXY_display2( Shift_map, chk_print, Detector, sample_para, search_Range_2D, d_Range_2D );
