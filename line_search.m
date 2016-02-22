function [ Al_line, Ni_line, Al_Abso, Ni_Abso, Absrp_line ] = line_search( chkXY, tot_Det_num, search_Deg, d_Deg, sample_para, holder_para, holder_frame_para, angle_search, SpuriousX)
%1D tilt series calculation
%Weizong Xu, March, 2015, wxu4@ncsu.edu

if (abs(chkXY-2)<0.0001)
    for i=1:tot_Det_num
        [Al_line(:,:,i), Ni_line(:,:,i)] = TiltY_search_parallel(search_Deg, d_Deg, sample_para, holder_para, holder_frame_para, angle_search(:,:,i), SpuriousX(:,i));
    end
else
    for i=1:tot_Det_num
        [Al_line(:,:,i), Ni_line(:,:,i)] = TiltX_search_parallel(search_Deg, d_Deg, sample_para, holder_para, holder_frame_para, angle_search(:,:,i), SpuriousX(:,i), chkXY);
    end
end

for i=1:tot_Det_num
    tot_xray(i)=sum(angle_search(:,5,i));
    Al_Abso(:,1,i)=Al_line(:,1,i);
    Al_Abso(:,2,i)=Al_line(:,2,i)/tot_xray(i);
    Ni_Abso(:,1,i)=Ni_line(:,1,i);
    Ni_Abso(:,2,i)=Ni_line(:,2,i)/tot_xray(i);

    Absrp_line(:,1,i)=Al_line(:,1,i);
    Absrp_line(:,2,i)=Al_line(:,2,i)./Ni_line(:,2,i);
end




end