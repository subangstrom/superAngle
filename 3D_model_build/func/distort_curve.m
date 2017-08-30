function [ p ] = distort_curve( p, Input_curve)
%Weizong Xu, August, 2017

Radius=Input_curve.Radius;
chk_shape=Input_curve.chk_shape; %1 convex -1 concavo
Curve_range=Input_curve.Curve_range*(max(p(1,:))-min(p(1,:)));
y_damp=Input_curve.y_damp;


%x_l=max(mx);
x_l=Curve_range(2)-Curve_range(1);
r_x=x_l/2+Curve_range(1);
r_z=-sqrt(Radius^2-(x_l/2)^2);
y_min=min(p(2,:)); %model edge y=0;
y_max=max(p(2,:));
for i=1:length(p(1,:))
    mx=p(1,i);
    if (mx>=Curve_range(1) && mx<=Curve_range(2))
        %mx=0:100;
        dz=chk_shape*(sqrt(Radius^2-(mx-r_x)^2)+r_z);
        %figure;plot(mx,dz)
        damp_ratio=1-y_damp*(p(2,i)-y_min)/(y_max-y_min);
        if (damp_ratio<0)
            damp_ratio=0;
        end
        dz=dz*damp_ratio;
        p(3,i)=p(3,i)+dz;
    end
end

end

