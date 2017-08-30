function [ model_slice, sample_para ] = slice_model( model_sample, model_ether, model_connect, sample_para, ether_start_p, chk_display, chk_direction )
%Weizong Xu, August, 2017
%do some calculation to slice model, later it
% clear;clc;close all;
% load ('model_connect_3sample_5ether.mat')
%%initialization
%tic
%ether_start_p=[0,0,0];
%point=[50,50,5];
%get tetrahedron number its is within ether
% tic
% [ ~, p_tetra_num1] = reg_point_tetrahedron( ether_start_p, model_ether, -1 );
% toc
t_sort=model_ether.t_sort;
t_sort_s=model_sample.t_sort;
list=1:length(t_sort);
list=[sample_para.p_tetra_num list];
p=model_ether.p;
p_s=model_sample.p;
nei_list_coplane=model_ether.nei_list_coplane;
nei_list_coplane_s=model_sample.nei_list_coplane;
tetra_volume=model_ether.tetra_volume;
% tic
[ p_tetra_num ] = reg_point_tetrahedron_fast2( list, ether_start_p, p, t_sort, tetra_volume );
% toc
if (p_tetra_num==0)
    %try low accuracy search, fix bug
    [ p_tetra_num ] = reg_point_tetrahedron_fast_acc( list, ether_start_p, p, t_sort, tetra_volume, 1e-6 );
end
if (p_tetra_num==0)
    %out of specimen
    model_slice.Model_Thickness=0;
    return;
else
    sample_para.p_tetra_num=p_tetra_num;
end

line_direction=[0,0,-1];
if (chk_direction==-1)
    line_direction=[0,0,1];
end
%%get two points that line intersect with tetrahedron in sample from ether start poing
mesh_num=p_tetra_num;
line_point=ether_start_p;
chk_hit=0;
n_forward=0;
sample_start_p=[];
sample_start_num=[];


%four_neighbour=model_ether.nei_list_coplane.four_neighbour;
while (chk_hit==0 && isnan(n_forward)==0)
%     tic
%     [ p_forward1, n_forward1, ~, ~, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_ether );
%     toc
    [ p_forward, n_forward, ~, ~, ~ ] = intercept_line_tetra_fast( line_point, line_direction, mesh_num, p, t_sort, nei_list_coplane );
    sample_list=model_connect.tetra_reg_ether{mesh_num};
    if (isempty(sample_list)==0)
        for ii=1:length(sample_list);
            %[ p_forward11, n_forward11, p_backward11, n_backward11, chk_out ] = intercept_line_tetra( line_point, line_direction, sample_list(ii), model_sample );
            [ p_forward1, n_forward1, p_backward1, n_backward1, chk_out ] = intercept_line_tetra_fast( line_point, line_direction, sample_list(ii), p_s, t_sort_s, nei_list_coplane_s );
            if (chk_out==1 && isnan(n_forward1)==0)
                sample_start_num=n_forward1;
                sample_start_p=p_forward1;
                chk_hit=1;
                break;
            end
            if (chk_out==1 && isnan(n_backward1)==0) %bug fix if n_forward is nan
                sample_start_num=n_backward1;
                sample_start_p=p_backward1;
                chk_hit=1;
                break;
            end
        end
    end
    mesh_num=n_forward;
    line_point=p_forward;
end
%end

%%get line start and line end of the sample from the sample starting point
%tic
[dist_path, p_intercept, sample_list_all] = get_intercept_model(sample_start_p,sample_start_num,line_direction,model_sample,0); %0 get intercept in both direction; >0 in forward direction <0 backward direction
%toc
if (chk_direction==-1)
    temp=p_intercept.p_start;
    p_intercept.p_start=p_intercept.p_end;
    p_intercept.p_end=temp;
end

%%multislice sample into N slices (N+1 points)
N=sample_para.Slice_t;
dz_min=max((model_sample.p(3,:))-min(model_sample.p(3,:)))/N/1.25;
dz_input=dist_path/N;
if (dz_input<dz_min && sample_para.t_chk==0)
    N=max(ceil(dist_path/dz_min),25); %min 25 points for intergration
end

if (sample_para.t_chk == 1 || sample_para.t_chk==10) %set thickness same as sample_para inputed.
    Thickness = sample_para.Thickness;
end

if (sample_para.t_chk==0) %use model defined shape
    Thickness = dist_path/model_sample.geo_real_ratio;
end

if (sample_para.t_chk~=0 && sample_para.t_chk~=1 && sample_para.t_chk~=10) %same location
    Thickness = sample_para.Thickness/cos(TiltX*pi/180)/cos(TiltY*pi/180);
end

model_slice.RealModel_ratio=Thickness/dist_path; %multiply it to get real distance
model_slice.dt_Real=Thickness/N;
model_slice.dt_Model=dist_path/N;
model_slice.N=N;
model_slice.Real_Thickness=Thickness;
model_slice.Model_Thickness=dist_path;
%toc

add_on_coor=[];
model_slice.n=[]; %bug fix if no model_slice.n is registered in the following search
for ii=1:N %topest point is moved from this list so that size of p is same as N
    model_slice.p{ii,1}=p_intercept.p_start+(p_intercept.p_end-p_intercept.p_start)*ii/N;
    add_on_coor=[add_on_coor model_slice.p{ii,1}'];
    for jj=1:length(sample_list_all)
        chk=point_in_tetra(model_slice.p{ii,1}, model_sample,sample_list_all(jj));
        if (chk==1)
            model_slice.n{ii,1}=sample_list_all(jj);
        end
    end
end

if (length(model_slice.n)~=N)
    %disp('Error! Sliced number is not equal to N. Try to fix it')
    model_slice.n=[];
    add_on_coor=[];
    for ii=1:N %topest point is moved from this list so that size of p is same as N
        model_slice.p{ii,1}=p_intercept.p_start+(p_intercept.p_end-p_intercept.p_start)*ii/N;
        add_on_coor=[add_on_coor model_slice.p{ii,1}'];
        for jj=1:length(sample_list_all)
            chk=point_in_tetra_acc(model_slice.p{ii,1}, model_sample,sample_list_all(jj),1e-6);
            if (chk==1)
                model_slice.n{ii,1}=sample_list_all(jj);
            end
        end
    end
    if (length(model_slice.n)==N)
        %disp('Issue fixed.');
    else
        %disp('Still have this issue. Slice model failed!');
        if (~isempty(model_slice.n))
            disp('Error! Sliced number is not equal to N. This will cause calculation failure!')
        end
    end    
end

%bug fix if start point is below sample
if ((isempty(sample_start_p) || isempty(sample_start_num)) && chk_direction==1)
   [ model_slice ] = slice_model( model_sample, model_ether, model_connect, sample_para,ether_start_p, chk_display, -1 ); 
   chk_display=0;
end

if ((isempty(sample_start_p) || isempty(sample_start_num)) && chk_direction==-1)
    %disp('Error! Out of specimen region!')
    model_slice.Model_Thickness=0;
    return;
end
%% look at mesh slice
if (chk_display==1)
    look_p=sample_list_all;
    add_on_coor=[add_on_coor sample_start_p'];
    look_model_tetra( model_sample, look_p, add_on_coor );

    line1=ether_start_p(1);
    line2=ether_start_p(2);
    line31=ether_start_p(3)-20-model_slice.Model_Thickness/2;
    line32=ether_start_p(3)+20+model_slice.Model_Thickness/2;
    figure;pdeplot3D(model_sample.p,model_sample.t);
    hold on;
    plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r')
    hold off;
    disp(['Model thickness (distpath):',num2str(model_slice.Model_Thickness)])
    if (model_sample.geo_real_ratio~=1)
        disp(['Model thickness offset by geo_real_ratio (distpath):',num2str(model_slice.Model_Thickness/model_sample.geo_real_ratio)])
    end
    disp(['Real thickness (thickness):',num2str(model_slice.Real_Thickness)])
% u_value=rand(length(p),3)*0.5;
% u_value=ones(length(p),3)*0;u_value(3,1)=1;
% figure;pdeplot3D(p,t,'colormapdata',u_value(:,1)); 
end

if (chk_display==2)
    disp(['Model thickness (distpath):',num2str(model_slice.Model_Thickness)])
    if (model_sample.geo_real_ratio~=1)
        disp(['Model thickness offset by geo_real_ratio (distpath):',num2str(model_slice.Model_Thickness/model_sample.geo_real_ratio)])
    end
    disp(['Real thickness (thickness):',num2str(model_slice.Real_Thickness)])
end

end