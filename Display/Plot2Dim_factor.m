function [ ] = Plot2Dim_factor(input, title_in, plotwidth, xaxis, yaxis, color_in, factor)
%Display in 2D for xy (n,2) format, n is total number
%Weizong Feb. 2015
%color_in = [0,0.5,1];
%plot(input(:,1),input(:,2)*factor,plotstyle, 'color',color_in);
plot(input(:,1),input(:,2)*factor,'LineWidth',plotwidth, 'color',color_in);
title(title_in); 
%text(0,0,0,'Sample'); 
xlabel(xaxis);
ylabel(yaxis);
%grid on;
box on;

end

