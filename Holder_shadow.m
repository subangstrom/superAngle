function [ H_angle_out ] = Holder_shadow(sample_para, holder_para, holder_frame_para, angle_search_in)
%Main code to calculate holder shadow and absorption effect
%Weizong Xu, April, 2015, wxu4@ncsu.edu

chk = holder_para(3);
if (chk>0)
    [angle_search_level1] = Holder_grid_shadow_up(sample_para, holder_para, angle_search_in);  
    [angle_search_level2] = Holder_frame_shadow_fast_up(sample_para, holder_para, holder_frame_para, angle_search_level1);
    %angle_search_level2 = angle_search_level1;
    [angle_search_level3] = Holder_carrier_shadow_up_2(sample_para, holder_para, angle_search_level2);
    H_angle_out= angle_search_level3;
    %detector_illustration( angle_search_level2 );
else
    H_angle_out=angle_search_in;

end

end

