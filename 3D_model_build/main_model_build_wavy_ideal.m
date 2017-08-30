%% Build dual specimen-vacuum model
% The output can be used for prediction
% Run SuperAngle_script_specimen_geometry_XXX.m further
% Weizong Xu, 2017, August

clear;clc;close all;
%load model data
temp=load ('model_either_210_210_15_4mesh.mat', 'model_mesh_nei');
model_ether=temp.model_mesh_nei;
temp=load ('model_sample_200_200_2_3mesh.mat', 'model_mesh_nei');
model_sample=temp.model_mesh_nei;
clear temp;
%% shift center to (000)
shift_v=-model_ether.start_point;
perturb_v=[0.01 -0.01 0]; %set a small perturb number to avoid co-point
shift_v=shift_v+perturb_v;
shift_v(1)=shift_v(1)-50;
shift_v(2)=shift_v(2)-50+100;
shift_v(3)=shift_v(3)-5.24-0.5;
model_ether = model_shift( model_ether, shift_v );
shift_v=-model_sample.start_point;
shift_v(2)=shift_v(2)+100;
shift_v=shift_v-perturb_v;
model_sample = model_shift( model_sample, shift_v );

%%rotate
%if needed, underconsruction

%% distort the model -curvature
Input_wave.myc=100; %arb unit
Input_wave.amp=0.10;
Input_wave.Curve_range=[-1,1];
Input_wave.y_damp=2;
% model_sample.p = distort_curve_incline( model_sample.p, Input_wave);
% [model_sample] = cal_volume( model_sample);%must recalculate volume

Input_wave.w_length=200; %arb unit
Input_wave.amp=4;
Input_wave.phase=0;
Input_wave.Curve_range=[-1,1];
%Input_wave.y_damp=1;
Input_wave.y_damp=0;
model_sample.p = distort_curve_sinwave( model_sample.p, Input_wave);
[model_sample] = cal_volume( model_sample);%must recalculate volume
%% model Dilation for real sample (unit nm)
scale_bar=1;%25;
model_sample = model_dilation( model_sample, [1 1 1]*scale_bar);
model_ether = model_dilation( model_ether, [1 1 1.2]*scale_bar); %1.3 if not damp

%% connection check
if ( max(model_sample.p(1,:))>max(model_ether.p(1,:)) || ...
     max(model_sample.p(2,:))>max(model_ether.p(2,:)) || ...    
     max(model_sample.p(3,:))>max(model_ether.p(3,:)) || ...
     min(model_sample.p(1,:))<min(model_ether.p(1,:)) || ...
     min(model_sample.p(2,:))<min(model_ether.p(2,:)) || ...
     min(model_sample.p(3,:))<min(model_ether.p(3,:)) )
    disp('Connection check failed! Model sample is out of ether, please check!')
    return;
else
    disp('Connection check passed!')
end

%% look at model
figure;pdeplot3D(model_sample.p,model_sample.t); 
figure;pdeplot3D(model_sample.p,model_sample.t); 
hold on;
%figure;
pdeplot3D(model_ether.p,model_ether.t,'FaceAlpha',0.5); 
hold off;

%% MOST IMPORTANT, calculate neighbour list between two model
%model_connect.tetra_reg_sample -> neighbour list (tetrahedron) of one sample tetrahedron wih respect to ether (ether tegrahdron number list)
%model_connect.tetra_reg_ether -> neighbour list (tetrahedron) of one ether tetrahedron with respect to sample (sample tegrahdron number list)
%[ model_connect ] = connect_model( model_sample, model_ether );
%% step 1
[ model_connect ] = connect_model1( model_sample, model_ether );
%% step 2
[ model_connect ] = connect_model2( model_sample, model_ether, model_connect );
%% visual check if connection is done correctly (sample)
% look_num=5; %have problem at 11
% look_p=model_connect.tetra_reg_sample{look_num};
% add_on_coor=[];
% look_two_model_tetra( model_ether, look_p, model_sample, look_num, add_on_coor );
% 
% %% visual check if connection is done correctly (ether)
% look_num=[3493]; %have problem at 11
% look_p=model_connect.tetra_reg_ether{look_num};
% add_on_coor=[];
% look_two_model_tetra( model_sample, look_p, model_ether, look_num, add_on_coor );
%% save data
model_ether = reduce_model( model_ether);
model_sample = reduce_model( model_sample);
save ('3D_model_sinwave.mat','model_ether','model_sample','model_connect');