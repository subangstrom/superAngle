function [ Volume ] = Vcal_tetrahedron( P1,P2,P3,P4 )
%Weizong Xu, August, 2017

P21 = P2-P1;%coords(T(:,2),:)-coords(T(:,1),:);
P31 = P3-P1;%coords(T(:,3),:)-coords(T(:,1),:);
P41 = P4-P1;%coords(T(:,4),:)-coords(T(:,1),:);
%Volume = 1/6*abs(dot(P21,cross(P31,P41,2),2));
%Volume = 1/6*abs(det([P21;P31;P41]));

% A11=P21(1);A12=P21(2);A13=P21(3);
% A21=P31(1);A22=P31(2);A23=P31(3);
% A31=P41(1);A32=P41(2);A33=P41(3);
%fastest way is do manually
Volume =1/6*abs(P21(1)*(P31(2)*P41(3)-P31(3)*P41(2))-P21(2)*(P31(1)*P41(3)-P31(3)*P41(1))+P21(3)*(P31(1)*P41(2)-P31(2)*P41(1)));


end

