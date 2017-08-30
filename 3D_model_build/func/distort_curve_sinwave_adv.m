function [ p ] = distort_curve_sinwave_adv( p, Input_wave)
%Weizong Xu, August, 2017
%sigmoid function S(t)=1/(1+e^{-t}) is also considered.

w_length=Input_wave.w_length; %arb unit
amp=Input_wave.amp;
phase=Input_wave.phase; %rad unit
Curve_rangeX=Input_wave.Curve_rangeX*(max(p(1,:))-min(p(1,:)));
Curve_rangeY=Input_wave.Curve_rangeY;
s_center=Input_wave.s_center;
s_spread=Input_wave.s_spread; %typical 2-5;
y_damp=Input_wave.y_damp;

y_min=min(p(2,:)); %model edge y=0;
y_max=max(p(2,:))*Curve_rangeY;


for i=1:length(p(1,:))
    mx=p(1,i);
    if (mx>=Curve_rangeX(1) && mx<=Curve_rangeX(2))
        %mx=0:100;
        dz=amp*sin(mx/w_length*pi*2+phase*pi/180);
        %figure;plot(mx,dz)
        damp_ratio=1-y_damp*(p(2,i)-y_min)/(y_max-y_min);
        if (damp_ratio<0)
            damp_ratio=0;
        end
        if (p(2,i)<y_max || y_damp==0)
            dz=dz*damp_ratio*(1-1/(1+exp(-(p(2,i)-s_center)/s_spread)));
            p(3,i)=p(3,i)+dz;
        end
    end
end

end

