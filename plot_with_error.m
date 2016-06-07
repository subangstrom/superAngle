function [ ] = plot_with_error( data_X, data_Y, d_err,color_code,shape_code, shape_size, range_x, range_y, chk)
%UNTITLED3 Summary of this function goes here
% chk=1 single figure, others add on figure.
if (chk==1)
    figure; %new figure
    set (gca, 'LineWidth',1.5);
end
line_width=2;
%xx = data_X(1):1:data_X(end);
%yy = spline(data_X,data_Y,xx);
%p1=plot(data_X, data_Y, 'o', xx, yy);

%p1=plot(data_X, data_Y, 'o');
%lsline;

p1=plot(data_X, data_Y,'-');
set (p1, 'LineWidth',line_width)
set (p1,'Color', color_code);
xlim([range_x(1) range_x(2)]);
ylim([range_y(1) range_y(2)]);
hold on;

          
color_code_err=color_code*0.9;
line_width_error=1.5;
tick_size=30;
%check matlab version
v_a = version('-release');
v_b = str2double(v_a(regexp(v_a,'\d')));
if (v_b>2014 ||  strcmp( v_a, '2014b')==1)
    p2=terrorbar(data_X,data_Y,d_err,tick_size);
else
    p2=errorbar(data_X,data_Y,d_err);
    errorbar_tick_mod(p2,tick_size);
    %set(p2, 'linestyle', 'none');
end
    set (p2, 'LineWidth',line_width_error);
    set (p2, 'color',color_code_err);



scatter(data_X,data_Y,shape_size,'filled',shape_code,...
              'MarkerEdgeColor',color_code*0.65,...
              'MarkerFaceColor',color_code, ...
              'LineWidth',1.5);
          

% if (chk==1)
%     hold off;
% end

end

