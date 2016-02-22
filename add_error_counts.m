function [ ] = add_error_counts(color_code, input, er_type, input2)
%Add error bar
%Weizong Feb. 2016

color_code=color_code*0.5;
line_width=1.5;
tick_size=30;

if (er_type==1) %error from counts
x=input(:,1);
y=input(:,2);
e=3*sqrt(y);
%p1=errorbar(x,y,e,'k');
p1=errorbar(x,y,e);
errorbar_tick(p1,tick_size);
set(p1, 'linestyle', 'none');
set (p1, 'LineWidth',line_width);
set (p1, 'color',color_code);
%set(p1, 'marker', '+');
end

if (er_type==2) %error from total counts
    temp=size(input);
    x=input(:,1,1);
    y=zeros(temp(1),1);
    e=zeros(temp(1),1);
    for i=1:temp(3)
        y_temp=input(:,2,i);    
        y=y+y_temp;
        e_temp=3*sqrt(y_temp);
        e=e+e_temp.^2; %z=ax+-by error_total=sqrt(a^2*ex^2+b^2*ey^2)
    end
    e=sqrt(e);
%p1=errorbar(x,y,e,'k');
p1=errorbar(x,y,e);
errorbar_tick(p1,tick_size);
set(p1, 'linestyle', 'none');
set (p1, 'LineWidth',line_width);
set (p1, 'color',color_code);
%set(p1, 'marker', '+');
end

if (er_type==3) %error from ratio
    x=input(:,1);
    y_A=input(:,2);
    y_B=input2(:,2);
    e_A=3*sqrt(y_A);
    e_B=3*sqrt(y_B);
    y=y_A./y_B;
    e=y.*sqrt((e_A./y_A).^2+(e_B./y_B).^2);  %error divided  
%p1=errorbar(x,y,e,'k');
p1=errorbar(x,y,e);
errorbar_tick(p1,tick_size);
set(p1, 'linestyle', 'none');
set (p1, 'LineWidth',line_width);
set (p1, 'color',color_code);
%set(p1, 'marker', '+');
end


if (er_type==4) %error from ratio
    temp=size(input);
    x=input(:,1,1);
    y_A=zeros(temp(1),1);
    y_B=zeros(temp(1),1);
    e_A=zeros(temp(1),1);
    e_B=zeros(temp(1),1);
    e=zeros(temp(1),1);
    for i=1:temp(3)
        y_temp=input(:,2,i);    
        y_A=y_A+y_temp;
        e_temp=3*sqrt(y_temp);
        e_A=e_A+e_temp.^2; %z=ax+-by error_total=sqrt(a^2*ex^2+b^2*ey^2)
        
        y_temp=input2(:,2,i);    
        y_B=y_B+y_temp;
        e_temp=3*sqrt(y_temp);
        e_A=e_A+e_temp.^2; %z=ax+-by error_total=sqrt(a^2*ex^2+b^2*ey^2)
    end
    e_A=sqrt(e_A);
    e_B=sqrt(e_B);
    y=y_A./y_B;
    e=y.*sqrt((e_A./y_A).^2+(e_B./y_B).^2);    
%p1=errorbar(x,y,e,'k');
p1=errorbar(x,y,e);
errorbar_tick(p1,tick_size);
set(p1, 'linestyle', 'none');
set (p1, 'LineWidth',line_width);
set (p1, 'color',color_code);
%set(p1, 'marker', '+');
end

end

