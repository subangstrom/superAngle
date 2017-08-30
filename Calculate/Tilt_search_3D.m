function [ Result ] = Tilt_search_3D( Detector, search_TiltX, search_TiltY, d_TiltX, d_TiltY, sample_para, holder_para, SpuriousX, model_sample_in, model_ether_in, model_connect, ether_start_p)
%1D tilt series calculation
%Weizong Xu, March, 2015, wxu4@ncsu.edu
%tic
tot_Det_num=Detector.tot_Det_num;
%angle_search=Detector.angle_search;
tot_num_x=(search_TiltX(2)-search_TiltX(1))/d_TiltX+1;
tot_num_y=(search_TiltY(2)-search_TiltY(1))/d_TiltY+1;
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

%TiltX=sample_para.TiltX;
%TiltY=sample_para.TiltY;
search_TiltX1=search_TiltX(1);
search_TiltY2=search_TiltY(2);
%p_tetra_num=[];
%ppp=ProgressBar(tot_num_x);
parfor i=1:tot_num_x
    %i
    p_tetra_num=[]; %remove if single core
    for j=1:tot_num_y
        %j
        Tilt_X=(i-1)*d_TiltX+search_TiltX1;
        Tilt_Y=search_TiltY2-(j-1)*d_TiltY;

        [ A_omega_out, B_omega_out, A_counts_out, B_counts_out, p_tetra_num  ] = single_spot_3D_parallel(Tilt_X, Tilt_Y, Detector, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p,p_tetra_num);
        A_map(j,i,:)=A_omega_out; %i,j must inverse for the display (x/y will reverse in this case)
        A_map_counts(j,i,:)=A_counts_out;
        B_map(j,i,:)=B_omega_out;
        B_map_counts(j,i,:)=B_counts_out;
    end
    %ppp.progress;
end
%ppp.stop;
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