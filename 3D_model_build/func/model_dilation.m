function [ model ] = model_dilation( model, dilation_v )
%Weizong Xu, August, 2017

model.start_point(1)=model.start_point(1)*dilation_v(1);
model.start_point(2)=model.start_point(2)*dilation_v(2);
model.start_point(3)=model.start_point(3)*dilation_v(3);

model.p(1,:)=model.p(1,:)*dilation_v(1);
model.p(2,:)=model.p(2,:)*dilation_v(2);
model.p(3,:)=model.p(3,:)*dilation_v(3);
[model] = cal_volume( model);%must recalculate volume
disp ('Done with shift dilation.')
end

