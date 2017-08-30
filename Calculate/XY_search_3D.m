function [ Result ] = XY_search_3D( Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p)
%For 2D map calculation
%Weizong Xu, March, 2015, wxu4@ncsu.edu

if (sample_para.cal_chk==1)
    uiwait(msgbox('Unable to calculate if composition is not known, please check specimen input file.'));
    disp('Unable to calculate if composition is not known, please check specimen input file.');
    Result.A_map=0;
    Result.B_map=0;
    Result.ratio_map=0;
    return;
end

tot_Det_num=Detector.tot_Det_num;
angle_search=Detector.angle_search;
disp ('Search along both X and Y-tilt directions')
for i=1:tot_Det_num
    disp(['Working on Detector #',num2str(i)]); %Display detector number
    [A_map(:,:,i), B_map(:,:,i)] = Tilt_search_parallel_3D(search_Deg, d_Deg, sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), model_sample, model_ether, model_connect, ether_start_p);
end

Al_count_all=A_map(:,:,1);
Ni_count_all=B_map(:,:,1);
for i=2:tot_Det_num
    Al_count_all=Al_count_all+A_map(:,:,i);
    Ni_count_all=Ni_count_all+B_map(:,:,i);
end
ratio_map=Al_count_all./Ni_count_all;

Result.A_map=A_map;
Result.B_map=B_map;
Result.ratio_map=ratio_map;

end

