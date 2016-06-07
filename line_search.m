function [ Result ] = line_search( chkXY, Detector, search_Deg, d_Deg, sample_para, holder_para, SpuriousX)
%1D tilt series calculation
%Weizong Xu, March, 2015, wxu4@ncsu.edu

tot_Det_num=Detector.tot_Det_num;
angle_search=Detector.angle_search;
if (abs(chkXY-2)<0.0001)
    disp('Search along Y-tilt direction')
    for i=1:tot_Det_num
        disp(['Working on Detector #',num2str(i)]);
        [Al_line(:,:,i), Ni_line(:,:,i)] = TiltY_search_parallel(search_Deg, d_Deg, sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i));
    end
else
    disp('Search along X-tilt direction')
    for i=1:tot_Det_num
        disp(['Working on Detector #',num2str(i)]);
        [Al_line(:,:,i), Ni_line(:,:,i)] = TiltX_search_parallel(search_Deg, d_Deg, sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), chkXY);
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


Result.A_line=Al_line;
Result.B_line=Ni_line;
Result.A_Abso=Al_Abso;
Result.B_Abso=Ni_Abso;
Result.Absrp_line=Absrp_line;
end