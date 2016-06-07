function [ Result ] = positionXY_search( Detector, search_Range, d_Range, sample_para, holder_para, SpuriousX)
%Wobble parameters to get the sensitivity of the count ratio
%Weizong Xu, March, 2015.

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
disp ('Search along both X and Y directions')
for i=1:tot_Det_num
    disp(['Working on Detector #',num2str(i)]); %Display detector number
    [A_map(:,:,i), B_map(:,:,i)] = Position_search_parallel(search_Range, d_Range, sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i));
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

