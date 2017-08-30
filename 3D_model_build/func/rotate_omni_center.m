function [ data_out ] = rotate_omni_center( axis, theta, center, data_in)
%Rotation around axis of theta degree towards center
%Code by Weizong 2016

data_tmp=data_in;
%shift relative to the center
data_tmp(:,1)=data_tmp(:,1)-center(1);
data_tmp(:,2)=data_tmp(:,2)-center(2);
data_tmp(:,3)=data_tmp(:,3)-center(3);

RotM=rotate_omni(axis,theta);

data_tmp2=data_tmp*RotM;

data_out=data_tmp2;

data_out(:,1)=data_tmp2(:,1)+center(1);
data_out(:,2)=data_tmp2(:,2)+center(2);
data_out(:,3)=data_tmp2(:,3)+center(3);

end

