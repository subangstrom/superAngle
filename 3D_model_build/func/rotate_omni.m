function [ Rot_Matrix ] = rotate_omni( axis, theta )
%Rotation matrix generator for omnipotent condition
%Input arbitrary axis and rotation angle, then done.
%Code by Weizong 2011 in FORTRAN, rewritten to Matlab m format.


c=axis/norm(axis);

matrix1=zeros(3,3);
matrix1(1,1)=1;
matrix1(2,2)=1;
matrix1(3,3)=1;

sintheta=sin(theta/180.0*pi);  %!In degrees
costheta=cos(theta/180.0*pi);


matrix2=zeros(3,3);
for i=1:3
 for j=1:3
    matrix2(i,j)=c(i)*c(j);
 end
end

matrix3=zeros(3,3);
matrix3(1,2)=-c(3);
matrix3(2,1)=c(3);

matrix3(1,3)=c(2);
matrix3(3,1)=-c(2);

matrix3(2,3)=-c(1);
matrix3(3,2)=c(1);

Rot_Matrix=costheta*matrix1+(1-costheta)*matrix2+sintheta*matrix3;

end

