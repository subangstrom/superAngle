function [  ] = PositionXY_display2_counts( chk_print, A_map, B_map, ratio_map, tot_Det_num, angle_search, sample_para, search_Range_2D, d_Range_2D )
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

if (sum(A_map(:))==0 && sum(B_map(:)) ==0)
    return;
end

[ convert_factor_A,convert_factor_B, tempA,tempB ] = absolute_scale_factor( sample_para);
A_map_counts=A_map*convert_factor_A;
B_map_counts=B_map*convert_factor_B;

A_map_all=sum(A_map_counts,3);
B_map_all=sum(B_map_counts,3);


sym_A = get_element_name(sample_para(13),sample_para(14));
sym_B = get_element_name(sample_para(15),sample_para(16));
image_range=-search_Range_2D:d_Range_2D:search_Range_2D;

for i=1:tot_Det_num
    figure;
    set(gca,'FontSize',15)
    imagesc (image_range, -image_range, A_map_counts(:,:,i)), axis image;
    title_name = ['Counts of ', sym_A, ' Xray arriving at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (mm)','FontSize',15)
    ylabel('Specimen Position Y (mm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    tot_solid_Angle=sum (angle_search(:,5,i)); %total solid angle
    max_counts_A=tot_solid_Angle*convert_factor_A;
    colormap (maptype);
    
    hold on;
    [C,h]=contour (image_range, -image_range, A_map_counts(:,:,i),'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
end

    figure;
    set(gca,'FontSize',15)
    imagesc (image_range, -image_range, B_map_counts(:,:,i)), axis image;
    title_name = ['Counts of ', sym_B, ' Xray arriving at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (mm)','FontSize',15)
    ylabel('Specimen Position Y (mm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    tot_solid_Angle=sum (angle_search(:,5,i)); %total solid angle
    max_counts_B=tot_solid_Angle*convert_factor_B;
    colormap (maptype);
    hold on;
    [C,h]=contour (image_range, -image_range, B_map_counts(:,:,i),'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;

    if (chk_print==1)
        print(['map_2D_shift_',title_name],'-dpng','-r300');
    end
end

    figure;
    set(gca,'FontSize',15)    
    imagesc (image_range, -image_range, A_map_all), axis image;
    title_name = ['Counts of ', sym_A, ' Xray arriving at all detectors'];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    temp=sum(angle_search,3);
    tot_solid_Angle=sum(temp(:,5)); %total solid angle
    max_counts_A=tot_solid_Angle*convert_factor_A;
    colormap (maptype);    
    hold on;
    [C,h]=contour (image_range, -image_range, A_map_all,'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
end

    figure;
    set(gca,'FontSize',15)    
    imagesc (image_range, -image_range, B_map_all), axis image;
    title_name = ['Counts of ', sym_B, ' Xray arriving at all detectors'];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    temp=sum(angle_search,3);
    tot_solid_Angle=sum(temp(:,5)); %total solid angle
    max_counts_B=tot_solid_Angle*convert_factor_B;
    colormap (maptype);    
    hold on;
    [C,h]=contour (image_range, -image_range, B_map_all,'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
end


for i=1:tot_Det_num
    figure;
    set(gca,'FontSize',15)
    imagesc (image_range, -image_range, A_map(:,:,i)./B_map(:,:,i)*convert_factor_A/convert_factor_B), axis image;
    title_name = ['Calculated ', sym_A, ' vs ', sym_B, ' at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (mm)','FontSize',15)
    ylabel('Specimen Position Y (mm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'FontSize',15);
    colormap (maptype);
    hold on;
    [C,h]=contour (image_range, -image_range, A_map(:,:,i)./B_map(:,:,i)*convert_factor_A/convert_factor_B,'k');
    clabel(C,h,'FontSize',15);
    hold off;
    if (chk_print==1)
        print(['map_2D_shift_',title_name],'-dpng','-r300');
    end
end

figure
set(gca,'FontSize',15)
imagesc (image_range, -image_range, ratio_map*convert_factor_A/convert_factor_B), axis image;
title_name =['Calculated ', sym_A, ' vs ', sym_B, ' of all detectors'];
title (title_name);
xlabel('Specimen Position X (mm)','FontSize',15)
ylabel('Specimen Position Y (mm)','FontSize',15)
axis xy;
h=colorbar;
set(h,'FontSize',15);
colormap (maptype);
hold on;
[C,h]=contour (image_range, -image_range, ratio_map*convert_factor_A/convert_factor_B,'k');
clabel(C,h,'FontSize',15);
hold off;
if (chk_print==1)
print(['map_2D_shift_',title_name],'-dpng','-r300');
end


end

