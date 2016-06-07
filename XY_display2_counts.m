function [  ] = XY_display2_counts( Tilt_map, chk_print, Detector, sample_para, search_Deg_2D, d_Deg_2D )
%Display 2D map for correction coefficient due to X-ray absorption and/or holder shadowing 
%Weizong Xu, wxu4@ncsu.edu, April 2015

%check matlab version, if equal or newer than 2014a, use parula as colormap
%older version use jet as colormap
v_a = version('-release');
v_b = str2double(v_a(regexp(v_a,'\d')));
if (v_b>=2014)
    maptype='parula';
else
    maptype='jet';
end

A_map=Tilt_map.A_map;
B_map=Tilt_map.B_map;
ratio_map=Tilt_map.ratio_map;
if (sum(A_map(:))==0 && sum(B_map(:)) ==0)
    return;
end

[ convert_factor_A,convert_factor_B, ~,~ ] = absolute_scale_factor( sample_para );
A_map_counts=A_map*convert_factor_A;
B_map_counts=B_map*convert_factor_B;

A_map_all=sum(A_map_counts,3);
B_map_all=sum(B_map_counts,3);

sym_A = get_element_name(sample_para.EleA_num,sample_para.EleA_shell);
sym_B = get_element_name(sample_para.EleB_num,sample_para.EleB_shell);
image_range=-search_Deg_2D:d_Deg_2D:search_Deg_2D;

tot_Det_num=Detector.tot_Det_num;
angle_search=Detector.angle_search;
for i=1:tot_Det_num
    figure;
    set(gca,'FontSize',15)
    
    imagesc (image_range, -image_range, A_map_counts(:,:,i)), axis image;
    title_name = ['Counts of ', sym_A, ' X-ray arriving at Detector ',num2str(i)];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    tot_solid_Angle=sum (angle_search(:,5,i)); %total solid angle
    max_counts_A=tot_solid_Angle*convert_factor_A;
    caxis([0,max_counts_A])
    colormap (maptype);
    
    hold on;
    [C,h]=contour (image_range, -image_range, A_map_counts(:,:,i),'k');
    clabel(C,h,'FontSize',20,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_',title_name],'-dpng','-r300');
end

    figure;
    set(gca,'FontSize',15)
    imagesc (image_range, -image_range, B_map_counts(:,:,i)), axis image;
    title_name = ['Counts of ', sym_B, ' X-ray arriving at Detector ',num2str(i)];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    tot_solid_Angle=sum (angle_search(:,5,i)); %total solid angle
    max_counts_B=tot_solid_Angle*convert_factor_B;
    caxis([0,max_counts_B])
    colormap (maptype);
    hold on;
    [C,h]=contour (image_range, -image_range, B_map_counts(:,:,i),'k');
    clabel(C,h,'FontSize',20,'LabelSpacing',600);
    hold off;

    if (chk_print==1)
        print(['map_2D_',title_name],'-dpng','-r300');
    end
end

    figure;
    set(gca,'FontSize',15)    
    imagesc (image_range, -image_range, A_map_all), axis image;
    title_name = ['Counts of ', sym_A, ' X-ray arriving at all detectors'];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    temp=sum(angle_search,3);
    tot_solid_Angle=sum(temp(:,5)); %total solid angle
    max_counts_A=tot_solid_Angle*convert_factor_A;
    colormap (maptype);
    
    hold on;
    [C,h]=contour (image_range, -image_range, A_map_all,'k');
    clabel(C,h,'FontSize',20,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_',title_name],'-dpng','-r300');
end

    figure;
    set(gca,'FontSize',15)    
    imagesc (image_range, -image_range, B_map_all), axis image;
    title_name = ['Counts of ', sym_B, ' X-ray arriving at all detectors'];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    temp=sum(angle_search,3);
    tot_solid_Angle=sum(temp(:,5)); %total solid angle
    max_counts_B=tot_solid_Angle*convert_factor_B;
    colormap (maptype);
    
    hold on;
    [C,h]=contour (image_range, -image_range, B_map_all,'k');
    clabel(C,h,'FontSize',20,'LabelSpacing',600);
    hold off;
if (chk_print==1)
    print(['map_2D_',title_name],'-dpng','-r300');
end

%figure
for i=1:tot_Det_num
    figure;
    set(gca,'FontSize',15)
    imagesc (image_range, -image_range, A_map(:,:,i)./B_map(:,:,i)*convert_factor_A/convert_factor_B), axis image;
    title_name = ['Calculated ', sym_A, ' vs ', sym_B, ' counts ratio at Detector ',num2str(i)];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    axis xy;
    h=colorbar;
    set(h,'fontsize',15);
    caxis([0,1*convert_factor_A/convert_factor_B])
    colormap (maptype);
    hold on;
    [C,h]=contour (image_range, -image_range, A_map(:,:,i)./B_map(:,:,i)*convert_factor_A/convert_factor_B,'k');
    clabel(C,h,'FontSize',20);
    hold off;
    if (chk_print==1)
        print(['map_2D_',title_name],'-dpng','-r300');
    end
end

figure
%imshow (ratio_map,[]);
set(gca,'FontSize',15)
imagesc (image_range, -image_range, ratio_map*convert_factor_A/convert_factor_B), axis image;
title_name =['Calculated ', sym_A, ' vs ', sym_B, ' counts ratio of all detectors'];
title (title_name);
xlabel('Tilt X (degree)','FontSize',15)
ylabel('Tilt Y (degree)','FontSize',15)
axis xy;
h=colorbar;
set(h,'fontsize',15);
colormap (maptype);
hold on;
[C,h]=contour (image_range, -image_range, ratio_map*convert_factor_A/convert_factor_B,'k');
clabel(C,h,'FontSize',20);
hold off;
if (chk_print==1)
    print(['map_2D_',title_name],'-dpng','-r300');
end


end

