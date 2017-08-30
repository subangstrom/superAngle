clear
%load('Area_searchXY_sample_wedge_sinwave_withgrid_tchk1.mat')
load('Area_searchXY_sample_wedge_sinwave_nogrid_tchk0.mat')
%load('Area_searchXY_sample_wedge_sinwave_grid_tchk0.mat')
%load('Area_searchXY_sample_wedge_sinwave_nogrid_tchk0_thetaE14.mat')
%load('Area_searchXY_sample_wedge_sinwave_nogrid_tchk0_thetaE22.mat')
%angle_range=10; %within deg % output 2D map 
angle_range=0:2:30; %within deg %output line plot of all tilt range
num_i=0;

% det_list{1}=[1,2];
% det_list{2}=[1,3];
% det_list{3}=[1,4];
% det_list{4}=[2,3];
% det_list{5}=[2,4];
% det_list{6}=[3,4];
% det_list{7}=[1,2,3];
% det_list{8}=[1,2,4];
% det_list{9}=[1,3,4];
% det_list{10}=[2,3,4];

det_list{1}=[1,2];
det_list{2}=[1,3];
% det_list{3}=[1,4];
% det_list{4}=[2,3];
% det_list{5}=[2,4];
% det_list{6}=[3,4];
det_list{3}=[1,2,3];
% det_list{8}=[1,2,4];
% det_list{4}=[1,3,4];
det_list{4}=[2,3,4];

plot_x=angle_range';
plot_y_sum=zeros(length(angle_range),1);
plot_y_invid=cell(Detector.tot_Det_num,1);
for i=1:Detector.tot_Det_num
  plot_y_invid{i,1}=zeros(length(angle_range),1);
end
plot_y_det_comb=cell(length(det_list),1);
for i=1:length(det_list)
  plot_y_det_comb{i,1}=zeros(length(angle_range),1);
end


Atomic_weight_A = get_element_weight(sample_para.EleA_num);
Atomic_weight_B = get_element_weight(sample_para.EleB_num);
ideal_atomic_ratio = sample_para.Atomic_ratio;
ideal_weight_ratio = ideal_atomic_ratio*Atomic_weight_A/Atomic_weight_B;
ideal_at_per = ideal_atomic_ratio/(1+ideal_atomic_ratio);
ideal_wt_per = ideal_weight_ratio/(1+ideal_weight_ratio);
convert_factor_A=sample_para.convert_factor_A;
convert_factor_B=sample_para.convert_factor_B;
k_AB_ideal = sample_para.k_AB_ideal;
atomic_ratio_AB = sample_para.Atomic_ratio;
k_AB_ideal_atomic = k_AB_ideal/Atomic_weight_A*Atomic_weight_B;
Int_ratio_factor= atomic_ratio_AB/k_AB_ideal_atomic;

for angle_r=angle_range%angle_range(1):angle_range(length(angle_range))

if exist('Tilt_map_std_all','var')
	chk_version=2;
else
    if exist('Tilt_map_std','var')
        chk_version=1;
        Omega_A_std_all=sum(Tilt_map_std.A_map,3);
        Count_A_std_all=Omega_A_std_all*convert_factor_A;
        Omega_B_std_all=sum(Tilt_map_std.B_map,3);
        Count_B_std_all=Omega_B_std_all*convert_factor_B;
        Ratio_AB_std_all_Omega=Omega_A_std_all./Omega_B_std_all;
        correction_factor = k_AB_ideal./Ratio_AB_std_all_Omega;
        Ratio_AB_std_all=Count_A_std_all./Count_B_std_all;
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
E_image_det_comb=cell(length(det_list));%zeros(ly,lx);


for i=1:t_num
    
        if (chk_version==2)
        	if (Tilt_map_std_all.t_chk==0)
                Tilt_map_std=Tilt_map_std_all.Tilt_map_std{i,1};
                Omega_A_std_all=sum(Tilt_map_std.A_map,3);
                %Count_A_std_all=Omega_A_std_all*convert_factor_A;
                Omega_B_std_all=sum(Tilt_map_std.B_map,3);
                %Count_B_std_all=Omega_B_std_all*convert_factor_B;
                Ratio_AB_std_all_Omega=Omega_A_std_all./Omega_B_std_all;
                correction_factor = k_AB_ideal./Ratio_AB_std_all_Omega;
                %Ratio_AB_std_all=Count_A_std_all./Count_B_std_all;

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
        %Count_A_all=Omega_A_all*convert_factor_A;
        Omega_B_all=sum(Tilt_map.B_map,3);
        %Count_B_all=Omega_B_all*convert_factor_B;
        Ratio_AB_all_Omega=Omega_A_all./Omega_B_all;
        comp_ratio_out_all=Ratio_AB_all_Omega*Int_ratio_factor.*correction_factor;
        wt_comp_all=comp_ratio_out_all./(1+comp_ratio_out_all);
        wt_comp_diff_all=(wt_comp_all-ideal_wt_per)/ideal_wt_per*100;
        
        %Ratio_AB_all=Count_A_all./Count_B_all;
        %Ratio_AB_all_diff=(Ratio_AB_all-Ratio_AB_std_all)./Ratio_AB_std_all*100;
        %Ratio_AB_all_diff_range=Ratio_AB_all_diff;
        %t_rangeY=size(Ratio_AB_all_diff,1);
        %t_rangeX=size(Ratio_AB_all_diff,2);
        t_rangeY=size(wt_comp_diff_all,1);
        t_rangeX=size(wt_comp_diff_all,2);
        for ti=1:t_rangeY
            for tj=1:t_rangeX
                TiltX=(tj-(t_rangeX-1)/2-1)*d_Deg_2D;
                TiltY=(ti-(t_rangeY-1)/2-1)*d_Deg_2D;
                Tilt=acosd(cosd(TiltX)*cosd(TiltY));
                if (Tilt>angle_r)
                    %Ratio_AB_all_diff_range(ti,tj)=0;
                    wt_comp_diff_all(ti,tj)=nan;
                end
            end
        end
%        Ratio_AB_all_diff_range=Ratio_AB_all_diff(1+s_range:t_range-s_range,1+s_range:t_range-s_range);
        %Error_max=max(max(abs(Ratio_AB_all_diff_range)));
        Error_max=max(max(abs(wt_comp_diff_all)));
        M_image(coor_y,coor_x)=1;
        E_image(coor_y,coor_x)=Error_max;
            
            for j=1:Detector.tot_Det_num
%                 tmp_ratio=Tilt_map.A_map(:,:,j)./Tilt_map.B_map(:,:,j)*convert_factor_A/convert_factor_B;
%                 tmp_ratio_std=Tilt_map_std.A_map(:,:,j)./Tilt_map_std.B_map(:,:,j)*convert_factor_A/convert_factor_B;
                tmp_ratio=Tilt_map.A_map(:,:,j)./Tilt_map.B_map(:,:,j);
                tmp_ratio_std=Tilt_map_std.A_map(:,:,j)./Tilt_map_std.B_map(:,:,j);
                correction_factor = k_AB_ideal./tmp_ratio_std;
                comp_ratio_out=tmp_ratio*Int_ratio_factor.*correction_factor;
                wt_comp=comp_ratio_out./(1+comp_ratio_out);
                wt_comp_diff=(wt_comp-ideal_wt_per)/ideal_wt_per*100;            
                
                
                %tmp_ratio_diff=(tmp_ratio-tmp_ratio_std)./tmp_ratio_std*100;
                %tmp_ratio_diff_range=tmp_ratio_diff;
                t_rangeY=size(wt_comp_diff,1);
                t_rangeX=size(wt_comp_diff,2);
                for ti=1:t_rangeY
                    for tj=1:t_rangeX
                        TiltX=(tj-(t_rangeX-1)/2-1)*d_Deg_2D;
                        TiltY=(ti-(t_rangeY-1)/2-1)*d_Deg_2D;
                        Tilt=acosd(cosd(TiltX)*cosd(TiltY));
                        if (Tilt>angle_r)
                            %tmp_ratio_diff_range(ti,tj)=0;
                            wt_comp_diff(ti,tj)=nan;
                        end
                    end
                end
                
                %Error_max=max(max(abs(tmp_ratio_diff)));
                %t_range=size(tmp_ratio_diff,1);
                %tmp_ratio_diff_range=tmp_ratio_diff(1+s_range:t_range-s_range,1+s_range:t_range-s_range);
                %Error_max=max(max(abs(tmp_ratio_diff_range)));
                Error_max=max(max(abs(wt_comp_diff)));
                E_image_indv{j}(coor_y,coor_x)=Error_max;
            end
%% combination of detectors (4 detector)
%% det 1 and 2
          for jk=1:length(det_list)  
            Tilt_map_A_sum=zeros(size(Tilt_map.A_map,1),size(Tilt_map.A_map,2));Tilt_map_A_sum_std=Tilt_map_A_sum;
            Tilt_map_B_sum=zeros(size(Tilt_map.B_map,1),size(Tilt_map.B_map,2));Tilt_map_B_sum_std=Tilt_map_B_sum;
            for j=1:length(det_list{jk})
                Tilt_map_A_sum=Tilt_map_A_sum+Tilt_map.A_map(:,:,det_list{jk}(j));
                Tilt_map_B_sum=Tilt_map_B_sum+Tilt_map.B_map(:,:,det_list{jk}(j));
                Tilt_map_A_sum_std=Tilt_map_A_sum_std+Tilt_map_std.A_map(:,:,det_list{jk}(j));
                Tilt_map_B_sum_std=Tilt_map_B_sum_std+Tilt_map_std.B_map(:,:,det_list{jk}(j));                
            end
                tmp_ratio=Tilt_map_A_sum./Tilt_map_B_sum;
                tmp_ratio_std=Tilt_map_A_sum_std./Tilt_map_B_sum_std;
                correction_factor = k_AB_ideal./tmp_ratio_std;
                comp_ratio_out=tmp_ratio*Int_ratio_factor.*correction_factor;
                wt_comp=comp_ratio_out./(1+comp_ratio_out);
                wt_comp_diff=(wt_comp-ideal_wt_per)/ideal_wt_per*100;            
                
                
                t_rangeY=size(wt_comp_diff,1);
                t_rangeX=size(wt_comp_diff,2);
                for ti=1:t_rangeY
                    for tj=1:t_rangeX
                        TiltX=(tj-(t_rangeX-1)/2-1)*d_Deg_2D;
                        TiltY=(ti-(t_rangeY-1)/2-1)*d_Deg_2D;
                        Tilt=acosd(cosd(TiltX)*cosd(TiltY));
                        if (Tilt>angle_r)
                            wt_comp_diff(ti,tj)=nan;
                        end
                    end
                end
                
                Error_max=max(max(abs(wt_comp_diff)));
                E_image_det_comb{jk}(coor_y,coor_x)=Error_max;
          end
            
            
            
            
            c_num=c_num+1;
end

num_i=num_i+1;
disp(['Angle = ',num2str(angle_r), ' out of ',num2str(angle_range(length(angle_range))),' (deg)'])
for iii=1:Detector.tot_Det_num
    plot_y_invid{iii,1}(num_i,1)=max(max(E_image_indv{iii}));
end
plot_y_sum(num_i,1)=max(max(E_image));

for iii=1:length(det_list)
    plot_y_det_comb{iii,1}(num_i,1)=max(max(E_image_det_comb{iii}));
end

end
%%

if (length(angle_range)==1)
    disp (['Currently, ',num2str(c_num/t_num*100),'% calculation is done'])
    figure;imagesc(M_image);axis image;title('Availability of data');
    figure;imagesc(E_image, [0 20]);axis image;colorbar;
    colormap (parula);
    title('Deviation (%) of k-value of detector sum')
    print('Search_error_XY_2Dmap_all','-depsc','-tiff');
    for j=1:Detector.tot_Det_num
        tmp_image=E_image_indv{j};
        figure;imagesc(tmp_image, [0 50]);axis image;colorbar;
        colormap (parula);
        title(['Deviation (%) of k-value at detector #',num2str(j)])
        %print(['Search_error_XY_2Dmap_Det',num2str(j)],'-depsc','-tiff');
    end
    for j=1:length(det_list)
        figure;imagesc(E_image_det_comb{j}, [0 20]);axis image;colorbar;
        title(['Detector #', num2str(det_list{j})]);
        colormap (parula);
        print(['Search_error_XY_2Dmap_Det',num2str(det_list{j})],'-depsc','-tiff');
    end
else
%plot line profile
    c_select=[241,196,90;190,29,44;0,144,189;137,198,55;126,47,142;162,20,47;77,190,238;0,0,0]/255;
    c_select2=[190,29,44;0,144,189;209,152,18;73,106,30;162,20,47;77,190,238;241,196,90;126,47,142;0,0,0]/255;
    color_all=c_select(5,:);
    
    figure

    for j=1:Detector.tot_Det_num
        plot_y=plot_y_invid{j};
        plot(plot_x,plot_y,'color',c_select(j,:),'LineWidth',1.5)
        hold on;
    end
    
    for j=1:length(det_list)
        plot_y=plot_y_det_comb{j};
        if (length(det_list{j})==2)
            plot(plot_x,plot_y,'--','LineWidth',2,'color',c_select2(j,:))
        end
        if (length(det_list{j})==3)
            plot(plot_x,plot_y,':','LineWidth',2,'color',c_select2(j,:))
        end
    end
    plot(plot_x,plot_y_sum,'color',color_all,'LineWidth',4)
    %set(gca,'YScale','log');
    legend('Det 1','Det 2','Det 3','Det 4','Det 1+2','Det 1+3','Det 1+2+3','Det 2+3+4','Sum','Location','northwest')
    xlabel('Tilt angle (deg)')
    ylabel('Error (%)')
    xlim([min(angle_range), max(angle_range)])
    ylim([0,75])
    %print('Search_error_XY_map','-depsc','-tiff');
    %print('Search_error_XY_map','-dpng','-r300');
end