function [ model ] = blunt_edge( model, amp, spread, wedge_angle, chk_display )
%Weizong Xu, August, 2017

model_save=model;
a=amp;
v=spread; %spread
max_y=10*v;
u=4*v; %offset
max_before=max(model.p(3,:));
max_after=max(model.p(2,:))*sind(1.5)+a/(1+exp(-(max_y-u)/v));


for i=1:length(model.p(2,:))
    if (model.p(3,i)>1e-5)
        if (model.p(2,i)<max_y)
            model.p(3,i)=(model.p(3,i)+a/(1+exp(-(model.p(2,i)-u)/v)))*(model.p(3,i)/(model.p(2,i)*sind(wedge_angle)));
        else
            model.p(3,i)=(model.p(3,i)+a/(1+exp(-(max_y-u)/v)))*(model.p(3,i)/(model.p(2,i)*sind(wedge_angle)));
        end
    end
end

for i=1:length(model.p(3,:))
    if (model.p(3,i)>0)
        model.p(3,i)=model.p(3,i)*max_before/max_after;
    end
end

[model] = cal_volume( model);%must recalculate volume

%%%Reference
% %% define a wedge function
% x=0:0.1:200;
% y=x*sind(1.5)+200*sind(1.5);
% figure;plot(x,y)
% xlim([0 200])
% ylim([0 200])
% %% define a Sigmoid function
% a=20;
% v=10; %spread
% x=0:0.1:10*v;
% u=5*v; %offset
% Sy=a./(1+exp(-(x-u)/v));
% figure;plot(x,Sy)
%% display the effect, wedge + Sigmoid function
if (chk_display==1)
    a=amp;
    v=spread; %spread
    max_x=10*v;
    %x=0:0.1:10*v;
    x=min(model.p(2,:)):0.1:max(model.p(2,:));
    y=zeros(length(x));
    u=4*v; %offset
    for i=1:length(x)
        if x(i)<=max_x
            y(i)=x(i)*sind(wedge_angle)+a/(1+exp(-(x(i)-u)/v));
        else
            y(i)=x(i)*sind(wedge_angle)+a/(1+exp(-(max_x-u)/v));
        end       
    end
    y=y*max_before/y(length(y));
    figure;plot(x,y)
    %ylim([min(x) max(x)]/3)

    figure;pdeplot3D(model_save.p,model_save.t);
    figure;pdeplot3D(model.p,model.t);
end


end

