function [ model ] = model_shift( model, shiftv )
%Weizong Xu, August, 2017

%slightly shift model a little bit to avoid pointing corner exactly in some cases
model.start_point(1)=model.start_point(1)+shiftv(1);
model.start_point(2)=model.start_point(2)+shiftv(2);
model.start_point(3)=model.start_point(3)+shiftv(3);

model.p(1,:)=model.p(1,:)+shiftv(1);
model.p(2,:)=model.p(2,:)+shiftv(2);
model.p(3,:)=model.p(3,:)+shiftv(3);

disp ('Done with shift model.')

end

