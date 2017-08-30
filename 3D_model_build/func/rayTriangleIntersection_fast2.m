function [flag, t] = rayTriangleIntersection_fast2 (o, d, p0, p1, p2) 
% Ray/triangle intersection using the algorithm proposed by Möller and Trumbore (1997). 
% 
% Input: 
% o : origin. 
% d : direction. 
% p0, p1, p2: vertices of the triangle. 
% Output: 
% flag: (0) Reject, (1) Intersect. 
% u,v: barycentric coordinates. 
% t: distance from the ray origin. 
% Author: 
% Originally written by Jesus Mena, edited by David Berman (dberm22@gmail.com)

%epsilon = 0.00001;

e1 = p1-p0; 
e2 = p2-p0; 

q = [d(2)*e2(3)-d(3)*e2(2), d(3)*e2(1)-d(1)*e2(3), d(1)*e2(2)-d(2)*e2(1)]; %cross product 
a = e1(1)*q(1)+e1(2)*q(2)+e1(3)*q(3); % determinant of the matrix M

if (a>-0.0000001 && a<0.0000001) 
%disp('the vector is parallel to the plane (the intersection is at infinity') 
flag=0; 
% u=0; 
% v=0; 
t=0; 
return; 
end

f = 1/a; 
s = o-p0; 
u = f*(s(1)*q(1)+s(2)*q(2)+s(3)*q(3)); 

if (u<0.0) 
% the intersection is outside of the triangle 
flag=0; 
% u=0; 
% v=0; 
t=0;
return; 
end

r = [s(2)*e1(3)-s(3)*e1(2), s(3)*e1(1)-s(1)*e1(3), s(1)*e1(2)-s(2)*e1(1)]; 
v = f*(d(1)*r(1)+d(2)*r(2)+d(3)*r(3)); 

if (v<0.0 || u+v>1.0) 
% the intersection is outside of the triangle 
flag=0; 
% u=0; 
% v=0; 
t=0;
return; 
end

t = f*(e2(1)*r(1)+e2(2)*r(2)+e2(3)*r(3)); % verified! 
flag = 1; 
return; 


end
