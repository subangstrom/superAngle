function [ ] = Plot2Dim(input, title_in, plotstyle, xaxis, yaxis)
%Display in 2D for xy (n,2) format, n is total number
%Weizong Feb. 2015

plot(input(:,1),input(:,2),plotstyle, 'color',[0,0.5,1]);  
title(title_in); 
%text(0,0,0,'Sample'); 
xlabel(xaxis);
ylabel(yaxis);
grid on;

end

