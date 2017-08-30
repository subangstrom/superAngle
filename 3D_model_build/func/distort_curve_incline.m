function [ p ] = distort_curve_incline( p, Input_wave)
%Weizong Xu, August, 2017

myc=Input_wave.myc; %arb unit
amp=Input_wave.amp; %center to begin incline
Curve_range=Input_wave.Curve_range*(max(p(1,:))-min(p(1,:)));
y_damp=Input_wave.y_damp;

y_min=min(p(2,:)); %model edge y=0;
y_max=max(p(2,:));

for i=1:length(p(1,:))
    mx=p(1,i);
    my=p(2,i);
    if (mx>=Curve_range(1) && mx<=Curve_range(2))
        %mx=0:100;
        if (my-myc)<0
            myl=my-myc;
        else
            myl=0;
        end
        dz=amp*myl;
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

