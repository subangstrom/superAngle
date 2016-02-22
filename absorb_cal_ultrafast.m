function [ tot_A1, tot_A2] = absorb_cal_ultrafast (u1,u2, angle_precal, dt, Thickness, S_n)
%Calculate average distance of x-ray path
%Weizong Xu, Feb. 2015

[i_negtive]=find(angle_precal(:,1)<0);
angle_precal(i_negtive,:)=[];
%        Xray_num(i,2)=temp_num(1);
%i

%angle_precal(:,1)=S_n(3)./angle_precal(:,1);
%tt=dt:dt:Thickness;
%tot_A1=sum(angle_precal(:,2).*sum(exp(-u1*angle_precal(:,1)*tt)')');
%tot_A2=sum(angle_precal(:,3).*sum(exp(-u2*angle_precal(:,1)*tt)')');

tot_A1=0;
tot_A2=0;         
for tt=dt:dt:Thickness
tot_A1=tot_A1+sum(angle_precal(:,2).*exp(-(u1*(S_n(3)*tt)./angle_precal(:,1))));
tot_A2=tot_A2+sum(angle_precal(:,3).*exp(-(u2*(S_n(3)*tt)./angle_precal(:,1))));
%speed up tip: matrix manupilation is much faster than for loop

% %for ii=1:temp_num(1)
% %    Xray_dist = (S_n(3)*tt)/angle_precal(ii,1);
% %    %Avoid X-ray fully absorb by sample with some large tilt angle  
% %    if (Xray_dist>0)
% %    tot_A1=tot_A1+angle_precal(ii,2)*exp(-(u1*Xray_dist));
% %    tot_A2=tot_A2+angle_precal(ii,3)*exp(-(u2*Xray_dist));
% %    end
%    %if (Xray_dist<=0)
%    %    ii
%    %end
%
%
% %end
                
% %               [tot_Al_tt, tot_Ni_tt] = absorb_cal_fast(sample_para, angle_precal, tt, S_n);
%                %tot_Al=tot_Al+tot_A1;
%               %tot_Ni=tot_Ni+tot_A2;
end

end

