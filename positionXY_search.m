function [ A_map, B_map, ratio_map ] = positionXY_search( tot_Det_num, search_Range, d_Range, sample_para, holder_para, holder_frame_para, angle_search, SpuriousX)
%Wobble parameters to get the sensitivity of the count ratio
%Weizong Xu, March, 2015.

if (sample_para(12)==1)
    uiwait(msgbox('Unable to calculate if composition is not known, please check specimen input file.'));
    disp('Unable to calculate if composition is not known, please check specimen input file.');
    A_map=0;
    B_map=0;
    ratio_map=0;
    return;
end

for i=1:tot_Det_num
    Detector_num = i %Display detector number
    [A_map(:,:,i), B_map(:,:,i)] = Position_search_parallel(search_Range, d_Range, sample_para, holder_para, holder_frame_para, angle_search(:,:,i), SpuriousX(:,i));
end

Al_count_all=A_map(:,:,1);
Ni_count_all=B_map(:,:,1);
for i=2:tot_Det_num
    Al_count_all=Al_count_all+A_map(:,:,i);
    Ni_count_all=Ni_count_all+B_map(:,:,i);
end
ratio_map=Al_count_all./Ni_count_all;




end

