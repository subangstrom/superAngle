function [  ] = sample_PosXY_display2_counts_Titan( Pos_map, Detector, sample_para, parameters)
%Display 2D map for correction coefficient due to X-ray absorption and/or holder shadowing
%Along Specimen shift in X and Y direction.
%Weizong Xu, wxu4@ncsu.edu, March 2015

%check matlab version, if equal or newer than 2014a, use parula as colormap
%older version use jet as colormap
v_a = version('-release');
v_b = str2double(v_a(regexp(v_a,'\d')));
if (v_b>=2014)
    maptype='parula';
else
    maptype='jet';
end

search_PosX=parameters.search_PosX;
d_PosX=parameters.d_PosX;
search_PosY=parameters.search_PosY;
d_PosY=parameters.d_PosY;
chk_print=parameters.chk_print;
rot_correction=parameters.rot_correction;
corr_probe=parameters.corr_probe;
corr_time=parameters.corr_time;
corr_geo=parameters.corr_geo;

disp_indv_ratioMap=parameters.disp_indv_ratioMap;
auto_scale_ratioMap=parameters.auto_scale_ratioMap;

disp_indv_countsMap=parameters.disp_indv_countsMap;
auto_scale_countsMap=parameters.auto_scale_countsMap;
if (auto_scale_countsMap==0 && disp_indv_countsMap~=0)
    background_level_countsMap=parameters.background_level_countsMap;
    disp('Background level set for counts is under construction...')
end

A_map=Pos_map.A_map;
B_map=Pos_map.B_map;
A_map_counts=Pos_map.A_map_counts;
B_map_counts=Pos_map.B_map_counts;
ratio_map=Pos_map.ratio_map;

A_map = imrotate(A_map,rot_correction,'bilinear','crop');
B_map = imrotate(B_map,rot_correction,'bilinear','crop');
A_map_counts = imrotate(A_map_counts,rot_correction,'bilinear','crop');
B_map_counts = imrotate(B_map_counts,rot_correction,'bilinear','crop');
ratio_map = imrotate(ratio_map,rot_correction,'bilinear','crop');

A_map_counts=A_map_counts*corr_probe*corr_time*corr_geo; %set as a standard condition
B_map_counts=B_map_counts*corr_probe*corr_time*corr_geo;

if (sum(A_map(:))==0 && sum(B_map(:)) ==0)
    return;
end

[ convert_factor_A,convert_factor_B, ~,~ ] = absolute_scale_factor( sample_para);
% A_map_counts=A_map*convert_factor_A;
% B_map_counts=B_map*convert_factor_B;

A_map_all=sum(A_map_counts,3);
B_map_all=sum(B_map_counts,3);


sym_A = get_element_name(sample_para.EleA_num,sample_para.EleA_shell);
sym_B = get_element_name(sample_para.EleB_num,sample_para.EleB_shell);
image_rangeX=search_PosX(1):d_PosX:search_PosX(2);
image_rangeY=search_PosY(2):-d_PosY:search_PosY(1);
tot_Det_num=Detector.tot_Det_num;
%angle_search=Detector.angle_search;

%%**********display counts map**************
if (disp_indv_countsMap~=0)
    
for i=1:tot_Det_num
    figure;
    set(gca,'FontSize',15)
    imagesc (image_rangeX, image_rangeY, A_map_counts(:,:,i)), axis image;
    title_name = ['Counts of ', sym_A, ' Xray arriving at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (mm)','FontSize',15)
    ylabel('Specimen Position Y (mm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    %tot_solid_Angle=sum (angle_search(:,5,i)); %total solid angle
    %max_counts_A=tot_solid_Angle*convert_factor_A;
    colormap (maptype);
    if (auto_scale_countsMap==0)
        caxis(parameters.Z_scale_indvA_countsMap);
    end
    hold on;
    [C,h]=contour (image_rangeX, image_rangeY, A_map_counts(:,:,i),'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
    print(['map_2D_shift_',title_name],'-depsc','-tiff');
end

    figure;
    set(gca,'FontSize',15)
    imagesc (image_rangeX, image_rangeY, B_map_counts(:,:,i)), axis image;
    title_name = ['Counts of ', sym_B, ' Xray arriving at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (nm)','FontSize',15)
    ylabel('Specimen Position Y (nm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    %tot_solid_Angle=sum (angle_search(:,5,i)); %total solid angle
    %max_counts_B=tot_solid_Angle*convert_factor_B;
    colormap (maptype);
    if (auto_scale_countsMap==0)
        caxis(parameters.Z_scale_indvB_countsMap);
    end
    hold on;
    [C,h]=contour (image_rangeX, image_rangeY, B_map_counts(:,:,i),'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;

    if (chk_print==1)
        print(['map_2D_shift_',title_name],'-dpng','-r300');
        print(['map_2D_shift_',title_name],'-depsc','-tiff');
    end
end

    figure;
    set(gca,'FontSize',15)    
    imagesc (image_rangeX, image_rangeY, A_map_all), axis image;
    title_name = ['Counts of ', sym_A, ' Xray arriving at all detectors'];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    %temp=sum(angle_search,3);
    %tot_solid_Angle=sum(temp(:,5)); %total solid angle
    %max_counts_A=tot_solid_Angle*convert_factor_A;
    colormap (maptype);
    if (auto_scale_countsMap==0)
        caxis(parameters.Z_scale_sumA_countsMap);
    end
    hold on;
    [C,h]=contour (image_rangeX, image_rangeY, A_map_all,'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
    print(['map_2D_shift_',title_name],'-depsc','-tiff');
end

    figure;
    set(gca,'FontSize',15)    
    imagesc (image_rangeX, image_rangeY, B_map_all), axis image;
    title_name = ['Counts of ', sym_B, ' Xray arriving at all detectors'];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    %temp=sum(angle_search,3);
    %tot_solid_Angle=sum(temp(:,5)); %total solid angle
    %max_counts_B=tot_solid_Angle*convert_factor_B;
    colormap (maptype);
    if (auto_scale_countsMap==0)
        caxis(parameters.Z_scale_sumB_countsMap);
    end
    hold on;
    [C,h]=contour (image_rangeX, image_rangeY, B_map_all,'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
    print(['map_2D_shift_',title_name],'-depsc','-tiff');
end

end
%%**********display ratio map**************
if (disp_indv_ratioMap~=0)
    
for i=1:tot_Det_num
    ratio_map_tmp=A_map(:,:,i)./B_map(:,:,i)*convert_factor_A/convert_factor_B;
    if (auto_scale_ratioMap==0)
        ratio_map_tmp(isnan(ratio_map_tmp))=parameters.Z_scale_ratioMap(1)+parameters.background_level_ratioMap*(parameters.Z_scale_ratioMap(2)-parameters.Z_scale_ratioMap(1));
    end
    figure;
    set(gca,'FontSize',15)
    imagesc (image_rangeX, image_rangeY, ratio_map_tmp), axis image;
    title_name = ['Calculated ', sym_A, ' vs ', sym_B, ' at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (nm)','FontSize',15)
    ylabel('Specimen Position Y (nm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    colormap (maptype);
    if (auto_scale_ratioMap==0)
        caxis(parameters.Z_scale_ratioMap);
    end
    hold on;
    [C,h]=contour (image_rangeX, image_rangeY, A_map(:,:,i)./B_map(:,:,i)*convert_factor_A/convert_factor_B,'k');
    clabel(C,h,'FontSize',15);
    hold off;
    if (chk_print==1)
        print(['map_2D_shift_',title_name],'-dpng','-r300');
        print(['map_2D_shift_',title_name],'-depsc','-tiff');
    end
end

figure
ratio_map_tmp=ratio_map;
if (auto_scale_ratioMap==0)
	ratio_map_tmp(isnan(ratio_map_tmp))=parameters.Z_scale_ratioMap(1)+parameters.background_level_ratioMap*(parameters.Z_scale_ratioMap(2)-parameters.Z_scale_ratioMap(1));
end
set(gca,'FontSize',15)
imagesc (image_rangeX, image_rangeY, ratio_map_tmp), axis image;
title_name =['Calculated ', sym_A, ' vs ', sym_B, ' of all detectors'];
title (title_name);
xlabel('Specimen Position X (nm)','FontSize',15)
ylabel('Specimen Position Y (nm)','FontSize',15)
axis xy;
h=colorbar;
set(h,'FontSize',15);
colormap (maptype);
if (auto_scale_ratioMap==0)
    caxis(parameters.Z_scale_ratioMap);
end
hold on;
[C,h]=contour (image_rangeX, image_rangeY, ratio_map,'k');
clabel(C,h,'FontSize',15);
hold off;
if (chk_print==1)
print(['map_2D_shift_',title_name],'-dpng','-r300');
print(['map_2D_shift_',title_name],'-depsc','-tiff');
end

end

end

