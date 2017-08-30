function [ratio_map_diff_all, comp_map_wt_diff_all, comp_map_at_diff_all] = diff_display2D( chk_print, A_map, B_map, ratio_map, A_map_std, B_map_std, ratio_map_std, ...
    tot_Det_num, search_Deg, sample_para_std, sample_para, search_Deg_2D, d_Deg_2D)
%Display difference duing wobbler test
%Weizong Xu, April, 2015

%check matlab version, if equal or newer than 2014a, use parula as colormap
%older version use jet as colormap
v_a = version('-release');
v_b = str2double(v_a(regexp(v_a,'\d')));
if (v_b>=2014)
    maptype='parula';
else
    maptype='jet';
end

%A_map_diff = (A_map-A_map_std)./A_map_std;
%B_map_diff = (B_map-B_map_std)./B_map_std;

ratio_map_diff_all = (ratio_map-ratio_map_std)./ratio_map_std;
%comp_map_diff_all = ratio_map_diff_all;

sym_A = get_element_name(sample_para(13),sample_para(14));
sym_B = get_element_name(sample_para(15),sample_para(16));
Atomic_weight_A = get_element_weight(sample_para(13));
Atomic_weight_B = get_element_weight(sample_para(15));


k_AB_ideal = sample_para_std(20); % ideal k ratio in for composition in wt%
atomic_ratio_AB = sample_para_std(17);

comp_at_std = atomic_ratio_AB/ (1+ atomic_ratio_AB);
comp_wt_std = (atomic_ratio_AB*Atomic_weight_A/Atomic_weight_B);
comp_wt_std = comp_wt_std /(1+comp_wt_std);


comp_map_at = (A_map./B_map)./(A_map_std./B_map_std)*atomic_ratio_AB;
comp_map_at = comp_map_at ./ (1+comp_map_at);
comp_map_at_all = (ratio_map)./(ratio_map_std)*atomic_ratio_AB;
comp_map_at_all = comp_map_at_all ./ (1+comp_map_at_all);
comp_map_at_diff = (comp_map_at - comp_at_std) ./comp_at_std;
comp_map_at_diff_all = (comp_map_at_all - comp_at_std) ./comp_at_std;


comp_map_wt = (A_map./B_map)./(A_map_std./B_map_std)*(atomic_ratio_AB*Atomic_weight_A/Atomic_weight_B);
comp_map_wt = comp_map_wt ./ (1+comp_map_wt);
comp_map_wt_all = (ratio_map)./(ratio_map_std)*(atomic_ratio_AB*Atomic_weight_A/Atomic_weight_B);
comp_map_wt_all = comp_map_wt_all ./ (1+comp_map_wt_all);
comp_map_wt_diff = (comp_map_wt - comp_wt_std) ./comp_wt_std;
comp_map_wt_diff_all = (comp_map_wt_all - comp_wt_std) ./comp_wt_std;


figure
subplot(1,2,1)
imagesc (abs(ratio_map_diff_all)), axis image, axis off;
title (['Calculated the difference of ', sym_A, '/', sym_B, ' counts ratio (absolute value) of all detectors']);
colorbar;
%caxis([0.75,1])
colormap (maptype);
hold on;
[C,h]=contour (abs(ratio_map_diff_all),'k');
clabel(C,h,'FontSize',15);
hold off;

%figure
subplot(1,2,2)
imagesc (ratio_map_diff_all), axis image, axis off;
title (['Calculated the difference of ', sym_A, '/', sym_B, ' counts ratio of all detectors']);
colorbar;
%caxis([0.75,1])
colormap (maptype);
hold on;
[C,h]=contour (ratio_map_diff_all,'k');
clabel(C,h,'FontSize',15);
hold off;


figure;
for i=1:tot_Det_num
    subplot(2,floor((tot_Det_num+1)/2),i)
    %imshow (A_map(:,:,i)./B_map(:,:,i),[]);
    imagesc (comp_map_wt(:,:,i)), axis image, axis off;
    title_name = ['Norminal ', sym_A,  ' composition (wt%) at Detector ',num2str(i)];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    h=colorbar;
    set(h,'fontsize',15);
    %colorbar;
    colormap (maptype);
    hold on;
    [C,h]=contour (comp_map_wt(:,:,i),'k');
    clabel(C,h,'FontSize',15);
    hold off;
end

figure
subplot(1,2,1)
imagesc (comp_map_wt_all), axis image, axis off;
title (['Norminal ', sym_A, ' compositon (wt%) of all detectors']);
colorbar;
%caxis([0.75,1])
colormap (maptype);
hold on;
[C,h]=contour (comp_map_wt_all,'k');
clabel(C,h,'FontSize',15);
hold off;

%figure
subplot(1,2,2)
imagesc (comp_map_at_all), axis image, axis off;
title (['Norminal ', sym_A, ' compositon (at%) of all detectors']);
colorbar;
%caxis([0.75,1])
colormap (maptype);
hold on;
[C,h]=contour (comp_map_at_all,'k');
clabel(C,h,'FontSize',15);
hold off;

figure;
for i=1:tot_Det_num
    subplot(2,floor((tot_Det_num+1)/2),i)
    %imshow (A_map(:,:,i)./B_map(:,:,i),[]);
    imagesc (abs(comp_map_wt_diff(:,:,i))), axis image, axis off;
    title_name = ['Variation of ', sym_A, ' composition (wt% absolute value) at Detector ',num2str(i)];
    title (title_name);
    xlabel('Tilt X (degree)','FontSize',15)
    ylabel('Tilt Y (degree)','FontSize',15)
    h=colorbar;
    set(h,'fontsize',15);
    %colorbar;
    colormap (maptype);
    hold on;
    [C,h]=contour (abs(comp_map_wt_diff(:,:,i)),'k');
    clabel(C,h,'FontSize',15);
    hold off;
end

figure
subplot(1,2,1)
imagesc (abs(comp_map_wt_diff_all)), axis image, axis off;
title (['Variation of ', sym_A, ' compositon (wt% absolute value) of all detectors']);
colorbar;
%caxis([0.75,1])
colormap (maptype);
hold on;
[C,h]=contour (abs(comp_map_wt_diff_all),'k');
clabel(C,h,'FontSize',15);
hold off;

%figure
subplot(1,2,2)
imagesc (abs(comp_map_at_diff_all)), axis image, axis off;
title (['Variation of ', sym_A, ' compositon (at% absolute value) of all detectors']);
colorbar;
%caxis([0.75,1])
colormap (maptype);
hold on;
[C,h]=contour (abs(comp_map_at_diff_all),'k');
clabel(C,h,'FontSize',15);
hold off;

if (chk_print>0)
image_range=-search_Deg_2D:d_Deg_2D:search_Deg_2D;
figure
set(gca,'FontSize',15)
imagesc (image_range,image_range, abs(comp_map_wt_diff_all)), axis image;
title (['Variation of ', sym_A, ' compositon (wt% absolute value) of all detectors']);
xlabel('Tilt X (degree)','FontSize',15)
ylabel('Tilt Y (degree)','FontSize',15)
h=colorbar;
set(h,'fontsize',15);
%caxis([0.75,1])
colormap (maptype);
hold on;
[C,h]=contour (image_range,image_range, abs(comp_map_wt_diff_all),'k');
clabel(C,h,'FontSize',15);
hold off;
if (chk_print==1)
print('Fig_out_wobb_kratio_2Dmap','-depsc','-tiff');
print('Fig_out_wobb_kratio_2Dmap','-dpng','-r300');
end

end

end

