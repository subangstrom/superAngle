clear
%load('Area_searchXY_sample_wedge_sinwave_withgrid_tchk1.mat')
load('Area_searchXY_sample_wedge_sinwave_nogrid_tchk0.mat')
s_range=8; %12 only zero deg; 10:+-5deg; 8:+-10 deg; 6:+-15deg; 4:+-20deg; 2:+-25deg; 0:+-30deg;   


if exist('Tilt_map_std_all','var')
	chk_version=2;
else
    if exist('Tilt_map_std','var')
        chk_version=1;
        Omega_A_std_all=sum(Tilt_map_std.A_map,3);
        %Omega_A_std_all(:,1)=Omega_A_std_all(:,1)/Detector.tot_Det_num;
        Omega_B_std_all=sum(Tilt_map_std.B_map,3);
        %Omega_B_std_all(:,1)=Omega_B_std_all(:,1)/Detector.tot_Det_num;
        %Ratio_AB_std_all(:,1)=Omega_A_std_all(:,1);
        Ratio_AB_std_all=Omega_A_std_all./Omega_B_std_all;
    else
        disp('Incorrect file format. Please check!')
        return;
    end
end

c_num=0;
t_num=size(control_table,1);
Mx=sort(unique(control_table(:,1)));
My=sort(unique(control_table(:,2)));
lx=length(Mx);
dx=Mx(2)-Mx(1);
ly=length(My);
dy=My(2)-My(1);
M_image=zeros(ly,lx);
E_image=zeros(ly,lx);
E_image_indv=cell(Detector.tot_Det_num,1);


for i=1:t_num
    
        if (chk_version==2)
        	if (Tilt_map_std_all.t_chk==0)
            	Tilt_map_std=Tilt_map_std_all.Tilt_map_std{i,1};
                Omega_A_std_all=sum(Tilt_map_std.A_map,3);
                Omega_B_std_all=sum(Tilt_map_std.B_map,3);
                Ratio_AB_std_all=Omega_A_std_all./Omega_B_std_all;
            else
            	Tilt_map_std=Tilt_map_std_all.Tilt_map_std{1,1};
            end 
        end
        
        sx=control_table(i,1);
        sy=control_table(i,2);
        coor_x=(sx-Mx(1))/dx+1;
        coor_y=(sy-My(1))/dy+1;
        Tilt_map=Tilt_search_Result_POS{i,1};
        Omega_A_all=sum(Tilt_map.A_map,3);
        Omega_B_all=sum(Tilt_map.B_map,3);
        Ratio_AB_all=Omega_A_all./Omega_B_all;
            
        Ratio_AB_all_diff=(Ratio_AB_all-Ratio_AB_std_all)./Ratio_AB_std_all*100;
        %Error_max=max(max(abs(Ratio_AB_all_diff)));
        t_range=size(Ratio_AB_all_diff,1);
        Ratio_AB_all_diff_range=Ratio_AB_all_diff(1+s_range:t_range-s_range,1+s_range:t_range-s_range);
        Error_max=max(max(abs(Ratio_AB_all_diff_range)));
        
        M_image(coor_y,coor_x)=1;
        E_image(coor_y,coor_x)=Error_max;
            
            for j=1:Detector.tot_Det_num
                tmp_ratio=Tilt_map.A_map(:,:,j)./Tilt_map.B_map(:,:,j);
                tmp_ratio_std=Tilt_map_std.A_map(:,:,j)./Tilt_map_std.B_map(:,:,j);
                tmp_ratio_diff=(tmp_ratio-tmp_ratio_std)./tmp_ratio_std*100;
                %Error_max=max(max(abs(tmp_ratio_diff)));
                t_range=size(tmp_ratio_diff,1);
                tmp_ratio_diff_range=tmp_ratio_diff(1+s_range:t_range-s_range,1+s_range:t_range-s_range);
                Error_max=max(max(abs(tmp_ratio_diff_range)));
                E_image_indv{j}(coor_y,coor_x)=Error_max;
            end

            c_num=c_num+1;
end
disp (['Currently, ',num2str(c_num/t_num*100),'% calculation is done'])
figure;imagesc(M_image);axis image;title('Availability of data');
figure;imagesc(E_image, [0 15]);axis image;colorbar;
title('Deviation (%) of k-value of detector sum')
%print('Error_map_det_sum','-depsc','-tiff');
%print('Error_map_det_sum','-dpng','-r300');
for j=1:Detector.tot_Det_num
    tmp_image=E_image_indv{j};
    figure;imagesc(tmp_image, [0 100]);axis image;colorbar;
    title(['Deviation (%) of k-value at detector #',num2str(j)])
    %print(['Error_map_det_',num2str(j)'],'-depsc','-tiff');
    %print(['Error_map_det_',num2str(j)'],'-dpng','-r300');
end
