function [ ether_start_p_out, ratio_chk ] = thickness_const_model( ether_start_p0, model_sample, model_ether, model_connect, sample_para_tmp, control )
%Weizong Xu, August, 2017

%tic

t_sort=model_ether.t_sort;
list=1:length(t_sort);
p=model_ether.p;

tetra_volume=model_ether.tetra_volume;
% tic
[ p_tetra_num ] = reg_point_tetrahedron_fast2( list, ether_start_p0, p, t_sort, tetra_volume );
% toc
if (p_tetra_num==0)
    %try low accuracy search, fix bug
    [ p_tetra_num ] = reg_point_tetrahedron_fast_acc( list, ether_start_p0, p, t_sort, tetra_volume, 1e-6 );
end
if (p_tetra_num==0)
    %out of specimen
    return;
else
    sample_para_tmp.p_tetra_num=p_tetra_num;
end

search_step=control.search_step; %along y direction in current model setting
pos_acc=control.pos_acc;
chk_display=control.chk_display;

[ ratio_chk ] = get_model_thickness( model_sample, model_ether, model_connect, sample_para_tmp, ether_start_p0, chk_display );
if (ratio_chk<1) %real/model
    %disp('ratio_chk<1')
    ether_start_p1=ether_start_p0;
    ether_start_p1(2)=ether_start_p0(2)-min(search_step,abs(ether_start_p0(2)*0.95));
    [ ratio_chk1 ] = get_model_thickness( model_sample, model_ether, model_connect, sample_para_tmp, ether_start_p1, chk_display );
end
if (ratio_chk>1) %real/model
    %disp('ratio_chk>1')
    ether_start_p1=ether_start_p0;
    ether_start_p1(2)=ether_start_p0(2)+min(search_step);
    [ ratio_chk1 ] = get_model_thickness( model_sample, model_ether, model_connect, sample_para_tmp, ether_start_p1, chk_display );
end

if ((ratio_chk-1)*(ratio_chk1-1)>0)
    disp('Error! Initial moving step is too small, cannot find good point!')
    return;
end

while (abs(ether_start_p0(2)-ether_start_p1(2))>pos_acc)
    ether_start_p_mid=0.5*(ether_start_p0+ether_start_p1);
    [ ratio_chk_mid ] = get_model_thickness( model_sample, model_ether, model_connect, sample_para_tmp, ether_start_p_mid, chk_display );
    if ((ratio_chk_mid-1)*(ratio_chk-1)<0)
        ether_start_p1=ether_start_p_mid;
    else
        ether_start_p0=ether_start_p_mid;
    end
end
ether_start_p_out=0.5*(ether_start_p0+ether_start_p1);
[ ratio_chk ] = get_model_thickness( model_sample, model_ether, model_connect, sample_para_tmp, ether_start_p_out, chk_display );
%toc


end

