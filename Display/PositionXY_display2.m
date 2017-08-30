function [  ] = PositionXY_display2( Shift_map, chk_print, Detector, sample_para, search_Range_2D, d_Range_2D )
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

A_map=Shift_map.A_map;
B_map=Shift_map.B_map;
ratio_map=Shift_map.ratio_map;
if (sum(A_map(:))==0 && sum(B_map(:)) ==0)
    return;
end

sym_A = get_element_name(sample_para.EleA_num,sample_para.EleA_shell);
sym_B = get_element_name(sample_para.EleB_num,sample_para.EleB_shell);
image_range=-search_Range_2D:d_Range_2D:search_Range_2D;
A_map_all=sum(A_map,3);
B_map_all=sum(B_map,3);

tot_Det_num=Detector.tot_Det_num;
%angle_search=Detector.angle_search;
for i=1:tot_Det_num
    figure;
    set(gca,'FontSize',15)
    num_xray_Det=1;
    imagesc (image_range, -image_range, A_map(:,:,i)/num_xray_Det), axis image;
    title_name = ['Effective solid angle of ', sym_A, ' X-ray arriving at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (mm)','FontSize',15)
    ylabel('Specimen Position Y (mm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    colormap (maptype);
    
    hold on;
    [C,h]=contour (image_range, -image_range, A_map(:,:,i)/num_xray_Det,'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
    print(['map_2D_shift_',title_name],'-depsc','-tiff');
end
    figure;
    set(gca,'FontSize',15)
    imagesc (image_range, -image_range, B_map(:,:,i)/num_xray_Det), axis image;
    title_name = ['Effective solid angle of ', sym_B, ' X-ray arriving at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (mm)','FontSize',15)
    ylabel('Specimen Position Y (mm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    colormap (maptype);
    hold on;
    [C,h]=contour (image_range, -image_range, B_map(:,:,i)/num_xray_Det,'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;

    if (chk_print==1)
        print(['map_2D_shift_',title_name],'-dpng','-r300');
        print(['map_2D_shift_',title_name],'-depsc','-tiff');
    end
end

    figure;
    set(gca,'FontSize',15)    
    imagesc (image_range, -image_range, A_map_all), axis image;
    title_name = ['Effective solid angle of ', sym_A, ' X-ray arriving at all detectors'];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    colormap (maptype);    
    hold on;
    [C,h]=contour (image_range, -image_range, A_map_all,'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
    print(['map_2D_shift_',title_name],'-depsc','-tiff');
end

    figure;
    set(gca,'FontSize',15)    
    imagesc (image_range, -image_range, B_map_all), axis image;
    title_name = ['Effective solid angle of ', sym_B, ' X-ray arriving at all detectors'];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    colormap (maptype);    
    hold on;
    [C,h]=contour (image_range, -image_range, B_map_all,'k');
    clabel(C,h,'FontSize',15,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_shift_',title_name],'-dpng','-r300');
    print(['map_2D_shift_',title_name],'-depsc','-tiff');
end


for i=1:tot_Det_num
    figure;
    set(gca,'FontSize',15)
    imagesc (image_range, -image_range, A_map(:,:,i)./B_map(:,:,i)), axis image;
    title_name = ['Calculated ', sym_A, ' vs ', sym_B, ' counts ratio at Detector ',num2str(i)];
    title (title_name);
    xlabel('Specimen Position X (mm)','FontSize',15)
    ylabel('Specimen Position Y (mm)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    colormap (maptype);
    hold on;
    [C,h]=contour (image_range, -image_range, A_map(:,:,i)./B_map(:,:,i),'k');
    clabel(C,h,'FontSize',15);
    hold off;
    if (chk_print==1)
        print(['map_2D_shift_',title_name],'-dpng','-r300');
        print(['map_2D_shift_',title_name],'-depsc','-tiff');
    end
end

figure
set(gca,'FontSize',15)
imagesc (image_range, -image_range, ratio_map), axis image;
title_name =['Calculated ', sym_A, ' vs ', sym_B, ' counts ratio of all detectors'];
title (title_name);
xlabel('Specimen Position X (mm)','FontSize',15)
ylabel('Specimen Position Y (mm)','FontSize',15)
axis xy;
h=colorbar;
set(h,'fontsize',15);
colormap (maptype);
hold on;
[C,h]=contour (image_range, -image_range, ratio_map,'k');
clabel(C,h,'FontSize',15,'LabelSpacing',300);
hold off;
if (chk_print==1)
print(['map_2D_shift_',title_name],'-dpng','-r300');
print(['map_2D_shift_',title_name],'-depsc','-tiff');
end

end

