function [ tot_A1, tot_A2] = absorb_cal_fast (sample_para, angle_precal, t, S_n)
%Calculate average distance of x-ray path
%Weizong Xu, Feb. 2015

% Note: I=I0*exp(-ul)=I0*exp[-(u/p)*pl]
% For Ni3Al system
% p=7.2 g/cm3
% (u/p)Al-K_Ni3Al = 4311.95 cm2/g (u/p)*p=31046 cm-1
% (u/p)Ni-K_Ni3Al = 61.012 cm2/g (u/p)*p=439.3 cm-1
u1 = sample_para(7)*1e-7;
u2 = sample_para(8)*1e-7;

%[S_n] = sample_normal(sample_para);

%probe_arriving --> [0, 0, t];
%S_coor = [0, 0, t];

%X_theta = pi/6;
%X_phi = pi/4;

tempSize=size(angle_precal);
%L_S = zeros(tempSize(1),3);
tot_A1=0;
tot_A2=0;


%avg_A1=mean(exp(-(u1*(S_n(3)*t)./angle_precal(:,1))));
%avg_A2=mean(exp(-(u2*(S_n(3)*t)./angle_precal(:,1))));

for i=1:tempSize(1)
%    %X_theta = angle_search(i,1);
%    %X_phi= angle_search(i,2);
%    %probe_arriving = [0, 0, t];
%    %Xray_dist = (S_n(3)*t)/ ...
%    %            (S_n(1)*sin(X_theta)*cos(X_phi)+S_n(2)*sin(X_theta)*sin(X_phi)+S_n(3)*cos(X_theta));
 
    Xray_dist = (S_n(3)*t)/angle_precal(i,1);
    %if (Xray_dist >350)

    %    [Xray_dist, t]
        
    %end


    %Avoid X-ray fully absorb by sample with some large tilt angle  
    if (Xray_dist>0)  
    %L_S(i,1) = Xray_dist*sin(X_theta)*cos(X_phi);
    %L_S(i,2) = Xray_dist*sin(X_theta)*sin(X_phi);
    %L_S(i,3) = Xray_dist*cos(X_theta);
    tot_A1=tot_A1+angle_precal(i,2)*exp(-(u1*Xray_dist));
    tot_A2=tot_A2+angle_precal(i,3)*exp(-(u2*Xray_dist));
    end
    
     %if (Xray_dist<0)
     %Xray_dist
     %angle_precal(i,1)
     %pause
     %end

end

%avg_A1 = tot_A1/double(tempSize(1));
%avg_A2 = tot_A2/double(tempSize(1));

%Display
%Scatter3Dim(L_S, 'sample')

end

