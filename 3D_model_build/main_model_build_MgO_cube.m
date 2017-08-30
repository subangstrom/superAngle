%% Build dual specimen-vacuum model
% The output can be used for prediction
% Run SuperAngle_script_specimen_geometry_XXX.m further
% Weizong Xu, 2017, August

clear;clc;close all;
%load model data
temp=load ('model_cube_110_20mesh', 'model_mesh_nei');
model_ether=temp.model_mesh_nei;
temp=load ('model_cube_200_200_200_12.5mesh.mat', 'model_mesh_nei');
model_sample=temp.model_mesh_nei;
clear temp;


%% rotate model
% axis=[0,-1,0];theta=90;center=[0,0,0];
% [ p ] = rotate_omni_center( axis, theta, center, model_ether.p');
% model_ether.p=p';
% [ p ] = rotate_omni_center( axis, theta, center, model_sample.p');
% model_sample.p=p';

%% shift center to (000)
shift_v=-model_ether.start_point;
perturb_v=[0.01 -0.01 0]; %set a small perturb number to avoid co-point
shift_v=shift_v+perturb_v;
shift_v(1)=shift_v(1);
shift_v(2)=shift_v(2);
shift_v(3)=shift_v(3);
model_ether = model_shift( model_ether, shift_v );
shift_v=-model_sample.start_point;
shift_v=shift_v-perturb_v;
model_sample = model_shift( model_sample, shift_v );

%%rotate
%if needed, underconsruction
% % distort the model -curvature
% Input_curve.Radius=150;
% Input_curve.chk_shape=-1; %1 convex -1 concavo
% Input_curve.Curve_range=[-0.33/2,0.33/2];
% Input_curve.y_damp=1;
% model_sample.p = distort_curve( model_sample.p, Input_curve);
% % shift_v=[0 0 0];
% % model_sample = model_shift( model_sample, shift_v );
% [model_sample] = cal_volume( model_sample);%must recalculate volume
%% model ether Dilation
dia_x=195/200;
dia_y=205/200;
dia_z=190/200;
model_sample = model_dilation( model_sample, [dia_x dia_y dia_z]);
model_ether = model_dilation( model_ether, [dia_x*2 dia_y*2 dia_z*2]);
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
% look_num=[343]; %have problem at 11
% look_p=model_connect.tetra_reg_ether{look_num};
% add_on_coor=[];
% look_two_model_tetra( model_sample, look_p, model_ether, look_num, add_on_coor );
%% save data
model_ether = reduce_model( model_ether);
model_sample = reduce_model( model_sample);
save ('model_connect_cubic_195_205_190.mat','model_ether','model_sample','model_connect');