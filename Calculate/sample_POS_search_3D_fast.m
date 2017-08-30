function [ Result ] = sample_POS_search_3D_fast( Detector, search_PosX, search_PosY, d_PosX, d_PosY, sample_para, holder_para, SpuriousX, model_sample_in, model_ether_in, model_connect)
%1D tilt series calculation
%Weizong Xu, March, 2015, wxu4@ncsu.edu
%tic
tot_Det_num=Detector.tot_Det_num;
%angle_search=Detector.angle_search;
tot_num_x_chk=(search_PosX(2)-search_PosX(1))/d_PosX+1;
tot_num_y_chk=(search_PosY(2)-search_PosY(1))/d_PosY+1;
x_range=search_PosX(1):d_PosX:search_PosX(2);
y_range=search_PosY(2):-d_PosX:search_PosY(1);
tot_num_x=length(x_range);
tot_num_y=length(y_range);
if (tot_num_x_chk ~= tot_num_x) %check if range is set fully
    disp('WARNING! Search X area does not cover full range set in search_PosX.')
    disp(['New range is set as [',num2str(x_range(1)),',',num2str(x_range(tot_num_x)),'].'])
else
    disp(['Search range along X is inputed as [',num2str(x_range(1)),',',num2str(x_range(tot_num_x)),'].'])
end
if (tot_num_y_chk ~= tot_num_y)
    disp('WARNING! Search Y area does not cover full range set in search_PosY.')
    disp(['New range is set as [',num2str(y_range(tot_num_y)),',',num2str(y_range(1)),'].'])
else
    disp(['Search range along Y is inputed as [',num2str(y_range(tot_num_y)),',',num2str(y_range(1)),'].'])
end

A_map=zeros(tot_num_y,tot_num_x,tot_Det_num);
B_map=zeros(tot_num_y,tot_num_x,tot_Det_num);
A_map_counts=zeros(tot_num_y,tot_num_x,tot_Det_num);
B_map_counts=zeros(tot_num_y,tot_num_x,tot_Det_num);

model_sample.p=model_sample_in.p;
model_sample.t_sort=model_sample_in.t_sort;
model_sample.tetra_volume=model_sample_in.tetra_volume;
model_sample.nei_list_coplane.four_neighbour=model_sample_in.nei_list_coplane.four_neighbour;
model_sample.geo_real_ratio=model_sample_in.geo_real_ratio;

model_ether.p=model_ether_in.p;
model_ether.t_sort=model_ether_in.t_sort;
model_ether.tetra_volume=model_ether_in.tetra_volume;
model_ether.nei_list_coplane.four_neighbour=model_ether_in.nei_list_coplane.four_neighbour;

TiltX=sample_para.TiltX;
TiltY=sample_para.TiltY;
search_PosX1=search_PosX(1);
search_PosY2=search_PosY(2);
ppp=ProgressBar(tot_num_x);
parfor i=1:tot_num_x
    %i
    p_tetra_num=[];
    
    for j=1:tot_num_y
        %[i,j]
        POS_X=(i-1)*d_PosX+search_PosX1;
        POS_Y=search_PosY2-(j-1)*d_PosY;
        ether_start_p=[POS_X, POS_Y, 0];
        
        [ A_omega_out, B_omega_out, A_counts_out, B_counts_out, p_tetra_num ] = single_spot_3D_parallel(TiltX, TiltY, Detector, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p, p_tetra_num);
        A_map(j,i,:)=A_omega_out; %i,j must inverse for the display (x/y will reverse in this case)
        A_map_counts(j,i,:)=A_counts_out;
        B_map(j,i,:)=B_omega_out;
        B_map_counts(j,i,:)=B_counts_out;
    end
    ppp.progress;
end
ppp.stop;
%toc

Al_count_all=A_map_counts(:,:,1);
Ni_count_all=B_map_counts(:,:,1);
for i=2:tot_Det_num
    Al_count_all=Al_count_all+A_map_counts(:,:,i);
    Ni_count_all=Ni_count_all+B_map_counts(:,:,i);
end
ratio_map=Al_count_all./Ni_count_all;

Result.A_map=A_map;
Result.B_map=B_map;

Result.A_map_counts=A_map_counts;
Result.B_map_counts=B_map_counts;
Result.ratio_map=ratio_map;
end