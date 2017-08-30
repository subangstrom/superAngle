function [ ether_start_p, chk_contact ] = tilt_compensate_startpoint( model_sample, model_ether, model_connect, sample_para, ether_start_p, choice )
%For ether_start_point, compensate tilt shift
%Weizong Xu, August, 2017
%tic
[ model_slice ] = slice_model( model_sample, model_ether, model_connect, sample_para, ether_start_p, 0, 1);
%toc
if (model_slice.Model_Thickness>0)
    chk_contact=1;
    p_top=model_slice.p{1};
    p_bottom=model_slice.p{model_slice.N};
    z_mid=(p_top(3)+p_bottom(3))/2;

    switch choice
        case 1
            ether_start_p(3)=z_mid;
        case 2
            ether_start_p(3)=p_top(3)+0.01;
        case 3
            ether_start_p(3)=p_bottom(3)+0.01;
        otherwise
            %ether_start_p=ether_start_p;
    end
else
    chk_contact=0;
end

end

