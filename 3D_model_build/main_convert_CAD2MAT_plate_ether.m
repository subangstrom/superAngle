%% use matlab function to generate mesh from AutoCAD input
% Raw input: stl output from AutoCAD or other software
% Allow reshape of the model in this script
% The matlab output can be further used to build dual specimen-vacuum model
% Please run main_model_build_XXX.m to build such dual model
% Weizong Xu, August 2017

clear; close all; clc
model = createpde(3);% Create a PDE Model
importGeometry(model,'plate_210_210_15.stl');% Construct the Geometry
%
% Plot the geometry and turn on face labels. You will need the face
% labels to define the boundary conditions.
figure
pdegplot(model,'FaceLabels','on');
view(30,30);
title('Bracket with Face Labels')

%Create a Mesh
% Create a mesh that uses 10-node tetrahedral elements with quadratic
% interpolation functions. This element type is significantly more accurate
% than the linear interpolation (four-node) elements, particularly in
% elasticity analyses that involve bending.
bracketThickness = 4; % Thickness of horizontal plate with hole, meters
hmax = bracketThickness; % Maximum element length for a moderately fine mesh
%generateMesh(model,'Hmax',hmax,'GeometricOrder','quadratic');
generateMesh(model,'Hmax',hmax,'GeometricOrder','linear');
figure
pdeplot3D(model);

%% grab mesh data from model *most important part
[p,~,t] = model.Mesh.meshToPet();
p_save=p;t_save=t;
p(1,:)=p_save(1,:)*1;
p(2,:)=p_save(2,:)*1;
p(3,:)=p_save(3,:)*1;
% find the corner that is close to zero point
p_test=p';
p_test(:,4)=(p_test(:,1).^2+p_test(:,2).^2+p_test(:,3).^2);
p_test(:,5)=1:length(p);
p_test2=sortrows(p_test,4);
zero_num=p_test2(1,5);
p(1,:)=p(1,:)-p_test2(1,1);
p(2,:)=p(2,:)-p_test2(1,2);
p(3,:)=p(3,:)-p_test2(1,3);

% rotate model
% axis=[0,-1,0];theta=90;center=[0,0,0];
% [ p ] = rotate_omni_center( axis, theta, center, p');
% p=p';

% distort the model -sin wave based.
% mx=1:100;
% dz=sin(mx/100*pi);
% figure;plot(mx,dz)
% p(3,:)=p(3,:)+10*sin(p(1,:)/100*pi);

% distort the model -curvature
Input_curve.Radius=150;
Input_curve.chk_shape=-1; %1 convex -1 concavo
Input_curve.Curve_range=[0.25,0.75];
Input_curve.y_damp=1;
%p = distort_curve( p, Input_curve);
%% display without model
%u_value=rand(length(p),3)*0.5;
u_value=ones(length(p),3)*0;u_value(3,1)=1;
figure;pdeplot3D(p,t);
hold on;
plot3([50,50],[50,50],[-30,30],'LineWidth',2)
hold off;
%figure;pdeplot3D(p,t,'colormapdata',u_value(:,1)); 

%% category of tetrahedrons (by faces) from t
%generate neighbour list (4 faces) of each tetrahedral 
%[p,~,t] = model.Mesh.meshToPet();
model_mesh.t=t;
model_mesh.p=p;
model_mesh.filename='model_ether.mat';
[model_mesh_nei]=cal_neighbour_list(model_mesh, 1); %do calculation if 1or2, if 2 in addition, save to model_sample.mat if zero load file instead

%% look at mesh
%look_p=1;
%look_p=[799,801,358,575,867];
look_p=model_mesh_nei.nei_list_copoint.seperate{1,1};
%look_p=1:868;
add_on_coor=[];
look_model_tetra( model_mesh_nei, look_p, add_on_coor );

%% input volume in the data structure
[ model_mesh_nei ] = cal_volume( model_mesh_nei);

%% get tetragonal number from a point
tic
point=[55,55,2.5];
[ model_mesh_nei] = reg_point_tetrahedron( point, model_mesh_nei, 1 );
toc
%% test: calculate intercept of a line from inner point of a tetrahedron
tic
line_direction=[0,0,1];

line_point=model_mesh_nei.start_point;
mesh_num=model_mesh_nei.start_tetra_num;
while (isnan(mesh_num)==0)
    [  p_forward, n_forward, p_backward, n_backward, chk_out ] = intercept_line_tetra_fast( line_point, line_direction, mesh_num, model_mesh_nei.p, model_mesh_nei.t_sort, model_mesh_nei.nei_list_coplane );
    line_point=p_forward;
    mesh_num=n_forward;
end
line_end=line_point;

line_point=model_mesh_nei.start_point;
mesh_num=model_mesh_nei.start_tetra_num;
while (isnan(mesh_num)==0)
    [  p_forward, n_forward, p_backward, n_backward, chk_out ] = intercept_line_tetra_fast( line_point, line_direction, mesh_num, model_mesh_nei.p, model_mesh_nei.t_sort, model_mesh_nei.nei_list_coplane );
    line_point=p_backward;
    mesh_num=n_backward;
end
line_start=line_point;

dist=norm(line_end-line_start);
toc

%% save file
save ('model_either_210_210_15_4mesh.mat', 'model_mesh_nei');