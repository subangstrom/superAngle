function [ ratio_chk, thickness ] = get_model_thickness( model_sample, model_ether, model_connect, sample_para_tmp, ether_start_p0, chk_display )
%Weizong Xu, August, 2017

[ ether_start_p, chk_contact ] = tilt_compensate_startpoint( model_sample, model_ether, model_connect, sample_para_tmp, ether_start_p0, 1 );

if (chk_contact==1) %have interception with sample
    [ model_sample_tmp, model_ether_tmp ] = tilt_model( model_sample, model_ether, sample_para_tmp );
    ether_start_p = sample_normal(ether_start_p,sample_para_tmp);
    [ model_slice ] = slice_model( model_sample_tmp, model_ether_tmp, model_connect, sample_para_tmp, ether_start_p, chk_display, 1);
else
    disp('Beam does not intercept with sample, please check!')
end
ratio_chk=model_slice.RealModel_ratio*model_sample.geo_real_ratio;
thickness=model_slice.Model_Thickness/model_sample.geo_real_ratio;
%toc

end

