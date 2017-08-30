function [ Result ] = line_search_3D_fast( chkXY, Detector, search_Deg, d_Deg, sample_para_in, holder_para, SpuriousX, model_sample_in, model_ether_in, model_connect, ether_start_p)
%1D tilt series calculation
%Weizong Xu, March, 2015, wxu4@ncsu.edu
%tic
tot_Det_num=Detector.tot_Det_num;
angle_search=Detector.angle_search;
tot_num=search_Deg/d_Deg*2+1;
Al_line=zeros(tot_num,2,tot_Det_num);
Ni_line=zeros(tot_num,2,tot_Det_num);

model_sample.p=model_sample_in.p;
model_sample.t_sort=model_sample_in.t_sort;
model_sample.tetra_volume=model_sample_in.tetra_volume;
model_sample.nei_list_coplane.four_neighbour=model_sample_in.nei_list_coplane.four_neighbour;
model_sample.geo_real_ratio=model_sample_in.geo_real_ratio;

model_ether.p=model_ether_in.p;
model_ether.t_sort=model_ether_in.t_sort;
model_ether.tetra_volume=model_ether_in.tetra_volume;
model_ether.nei_list_coplane.four_neighbour=model_ether_in.nei_list_coplane.four_neighbour;


for i=1:tot_num
    Al_line(i,1,:)=(i-1)*d_Deg-search_Deg;
    Ni_line(i,1,:)=(i-1)*d_Deg-search_Deg;
end
ppp=ProgressBar(tot_num);
parfor i=1:tot_num
    %i
    sample_para=sample_para_in;
    p_tetra_num=sample_para.p_tetra_num;
    if (chkXY==2)    
        TiltY=(i-1)*d_Deg-search_Deg;
        TiltX=sample_para.TiltX;
    else
        TiltX=(i-1)*d_Deg-search_Deg;
        TiltY=sample_para.TiltY;
    end
    %[Al_line(:,:,i), Ni_line(:,:,i)] = TiltY_search_parallel_3D_2(search_Deg, d_Deg, sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), ether_start_p, model_sample, model_ether, model_connect);
    [ A_omega_out, B_omega_out ] = single_spot_3D_parallel(TiltX, TiltY, Detector, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p, p_tetra_num);
    %Al_line(i,1,:)=(i-1)*d_Deg-search_Deg;
    Al_line(i,2,:)=A_omega_out;
    %Ni_line(i,1,:)=(i-1)*d_Deg-search_Deg;
    Ni_line(i,2,:)=B_omega_out;
    ppp.progress;
end
ppp.stop;
%toc


for i=1:tot_Det_num
    tot_xray(i)=sum(angle_search(:,5,i));
    Al_Abso(:,1,i)=Al_line(:,1,i);
    Al_Abso(:,2,i)=Al_line(:,2,i)/tot_xray(i);
    Ni_Abso(:,1,i)=Ni_line(:,1,i);
    Ni_Abso(:,2,i)=Ni_line(:,2,i)/tot_xray(i);

    Absrp_line(:,1,i)=Al_line(:,1,i);
    Absrp_line(:,2,i)=Al_line(:,2,i)./Ni_line(:,2,i);
end


Result.A_line=Al_line;
Result.B_line=Ni_line;
Result.A_Abso=Al_Abso;
Result.B_Abso=Ni_Abso;
Result.Absrp_line=Absrp_line;
Result.location_coor=ether_start_p;
end