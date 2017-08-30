%% Basic initiate
clear;close all;clc;
%delete progress*.txt
load ('model_connect_200_wedge_sinwave_incline_damp.mat');
%% prepare for parallel running
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool(8);%('local',8);
end
disp(['Note: CPU cores for this run is ',num2str(poolobj.NumWorkers),'.'])
%% **********Setup Detector***************
dAngle = 0.5; %Accuracy for angular intergration, 0.5 is fairly good for speed, 0.2 is better for accuracy
[ Detector ] = Detector_input( 'detector.xlsx', dAngle );
% %*****Spurious Xray calculation ***********
%Be cautious!!! Set to zero if not sure
[ SpuriousX ] = SpuriousXray_setup( Detector.tot_Det_num, 0, 0 );
%% **********Setup specimen************
t_chk = 0; % 1 --> constant thickness t input from sample parameter, 0- use model defined thickness(const spot) (nm), 10- match thickness with model defined thickness (search spot), other value --> constant spot during tilt, i.e. t will change
[ sample_para ] = Specimen_setup( 'specimen_wedge_ideal.xlsx', t_chk );
%[ sample_para ] = Specimen_setup( 'specimen_wedge_real.xlsx', t_chk );
sample_para.Thickness=50; %forcely set specimen thickness as 50 nm
sample_para.Z_rotation=0; %degree
[ model_sample, model_ether ]= tilt_model_precal( 0, 0, sample_para.Z_rotation, model_sample, model_ether);
model_sample.geo_real_ratio=1/25; % thickness calculation for nanoparticles, i.e. set sample size. %nm in unit
sample_para.chk_2nd=1; %1 consider 2nd interception (25-100x slow), other not consider
%**********Setup holder**************
chk_Shadow = 2;
%>0 and <>2) consider simple holder shadowing (circular contours and fully cutoff)
%2) consider holder shadowing by FEI low background double tilt HIVis holder
%<=0 value) No shadowing effect
[ holder_para ] = Holder_setup( chk_Shadow, 'holder_FEI_LB.xlsx', sample_para );
%[ holder_para ] = Holder_setup( chk_Shadow, 'holder_FEI_LB_demo_sample.xlsx', sample_para );
exp_file = 'exp_Ni3Al_wedge.xlsx';
chk_print = 2; %1 output high quality image to file, others>0 display but not output, <0 not display

search_Deg=30;
d_Deg = 1;
chkXY = 1; %2-Search Y, others search X
%% set tiltXY parameter
search_Deg_2D=30;
d_Deg_2D = 2.5;
search_TiltX=[-search_Deg_2D,search_Deg_2D];
search_TiltY=[-search_Deg_2D,search_Deg_2D];
d_TiltX=d_Deg_2D;
d_TiltY=d_Deg_2D;
chk_print = 2;
sample_para.Slice_t=25;

%% set parameter ideal in come cases
sample_para.Deviation_angle_X=0;
sample_para.Deviation_angle_Y=0;
sample_para.TiltX=0;
sample_para.TiltY=0;
sample_para.POSX=0;
sample_para.POSY=0;
%% Save a standard data from ideal plate shape (TiltX only)
% [ line_search_Result_std ] = line_search( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX);
% save('line_search_Result_std.mat','line_search_Result_std', 'chkXY', 'Detector', 'search_Deg', 'd_Deg', 'sample_para', 'holder_para', 'SpuriousX')
% line_display_Counts( line_search_Result_std, chk_print, chkXY, Detector, exp_file, search_Deg, sample_para );
%% Save a standard data from ideal plate shape (TiltXY)
%has some problem when t_chk=0; need fix
if (exist('Tilt_search_Result_std.mat','file'))
    disp('Skipped std calculation. Exist Tilt_search_Result_std.mat')
else
    disp('Run a standard tilt based on ideal plate shape')
    if (sample_para.t_chk~=0)
        [ Tilt_map_std ] = XY_search( Detector, search_Deg_2D, d_Deg_2D, sample_para, holder_para, SpuriousX);
        Tilt_map_std_all.Tilt_map_std{1,1}=Tilt_map_std;
        Tilt_map_std_all.t_chk=sample_para.t_chk;
        save('Tilt_search_Result_std.mat','Tilt_map_std', 'Detector', 'search_Deg_2D', 'd_Deg_2D', 'sample_para', 'holder_para', 'SpuriousX')
    else
        load('control.mat')
        for i=1:size(control_table,1)  %21 points
            sx=control_table(i,1);
            sy=control_table(i,2);
            ether_start_p0=[sx,sy,0]; %starting point of the model relative to ether (if sample shape varies, this could be useful)
            sample_para_tmp=sample_para;
            [ ~, thickness ] = get_model_thickness( model_sample, model_ether, model_connect, sample_para_tmp, ether_start_p0, 0 );
            sample_para_tmp.Thickness=thickness;
            disp (['#',num2str(i),' thickness=',num2str(thickness), 'nm'])
            [ Tilt_map_std ] = XY_search( Detector, search_Deg_2D, d_Deg_2D, sample_para_tmp, holder_para, SpuriousX);
            Tilt_map_std_all.Tilt_map_std{i,1}=Tilt_map_std;
            Tilt_map_std_all.Thickness{i,1}=thickness;
            Tilt_map_std_all.ether_start_p{i,1}=ether_start_p0;
        end
        Tilt_map_std_all.t_chk=0;
        save('Tilt_search_Result_std.mat','Tilt_map_std_all', 'Detector', 'search_Deg_2D', 'd_Deg_2D', 'sample_para', 'holder_para', 'SpuriousX')
    end
    %%XY_display2_counts( Tilt_map_std, chk_print, Detector, sample_para, search_Deg_2D, d_Deg_2D );
end
%*************************************************************
%% pickup one point and calculate
tot_cal=0;
rng('shuffle')
chk_p=0;chk_end=0;
while(chk_end==0)
load('control.mat', 'control_table')
 while(chk_p==0)
    rand_num=ceil(rand*size(control_table,1));
    if (control_table(rand_num,3)==1)
        sx=control_table(rand_num,1);
        sy=control_table(rand_num,2);
        %filename=['line_search_Result_',num2str(sx),'_',num2str(sy),'.mat'];
        filename=['Tilt_search_Result_',num2str(sx),'_',num2str(sy),'.mat'];
        if (exist(filename,'file'))
            control_table(rand_num,3)=0;
        else
            i=rand_num;
            chk_p=1;
            control_table(rand_num,3)=2; %pick up one node to calculate
            save ('control.mat','control_table');
        end
    end
    chk_end=1;
    for ii=1:size(control_table,1)
        if control_table(ii,3)==1
            chk_end=0;
        end
    end
    if (chk_p==0 && chk_end==1)
        disp('No more info needs to be processed. All set.');
        return; 
    end
 end
%for i=2:size(Data_cell,2)
    tic
    disp(['#',num2str(tot_cal+1)])
    sx=control_table(rand_num,1);
    sy=control_table(rand_num,2);
    ether_start_p=[sx,sy,0];
    %num_t=num_t+1;
    %disp(['#',num2str(num_t),' of ', num2str(((X_range(2)-X_range(1))/X_step+1)*((Y_range(2)-Y_range(1))/Y_step+1))])
    disp(num2str(ether_start_p));
    %[ line_search_Result ] = line_search_3D_fast( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
    [ Tilt_map ] = Tilt_search_3D( Detector, search_TiltX, search_TiltY, d_TiltX, d_TiltY, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
    tot_cal=tot_cal+1;
    %filename=['line_search_Result_',num2str(sx),'_',num2str(sy),'.mat'];
    %save(filename,'line_search_Result','ether_start_p')
    filename=['Tilt_search_Result_',num2str(sx),'_',num2str(sy),'.mat'];
    save(filename,'Tilt_map','ether_start_p')
    disp(['Search for position x=',num2str(sx),' y=', num2str(sy), ' is done.']);
    
    load ('control.mat','control_table');
    control_table(rand_num,3)=0;
    save ('control.mat','control_table');
    chk_p=0;
    toc
    if (tot_cal>=50)
        disp('Max number of work (50) for this node has been reached. Shutdown')
        %delete(poolobj);
        return;
    end
end
disp('No more info needs to be processed. All set.');
%delete(poolobj);

%% test and display
% sx=600;
% sy=2300;
% ether_start_p=[sx,sy,0];
% disp(['Working on position x=',num2str(sx),' y=', num2str(sy),'...']);
% tic
% [ Tilt_map ] = Tilt_search_3D( Detector, search_TiltX, search_TiltY, d_TiltX, d_TiltY, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p);
% toc
%%
 XY_display2_counts_3D( Tilt_map, 2, Detector, sample_para, search_Deg_2D, d_Deg_2D );
