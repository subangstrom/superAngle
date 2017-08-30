clear
load('control.mat','control_table')
Mx=sort(unique(control_table(:,1)));
My=sort(unique(control_table(:,2)));
lx=length(Mx);
dx=Mx(2)-Mx(1);
ly=length(My);
dy=My(2)-My(1);
M_image=zeros(ly,lx);
%check control table
n_node=0;
n_work=0;
t_num=size(control_table,1);
for i=1:t_num
    if control_table(i,3)==2
        n_node=n_node+1;
    end
    
    if control_table(i,3)==0
        n_work=n_work+1;
    end    
end

%check processed file
n_file=0;
for i=1:t_num
        sx=control_table(i,1);
        sy=control_table(i,2);
        coor_x=(sx-Mx(1))/dx+1;
        coor_y=(sy-My(1))/dy+1;
        filename=['Tilt_search_Result_',num2str(sx),'_',num2str(sy),'.mat'];
        if exist(filename, 'file')
            n_file=n_file+1;
            M_image(coor_y,coor_x)=1;
        end
end

if (n_file==n_work)
    disp ('control.mat file looks good!.')
else
    disp ('Control table record is not correct, regarding #points has been done. Please run repair file')
end


disp (['Currently, ',num2str(n_file/t_num*100),'% calculation is done'])
figure;imagesc(M_image);axis image;title('Availability of data');colormap(jet);
if (n_node>=1)
    disp (['Currently, ',num2str(n_node),' nodes are working on this calculation.'])
else
    if (n_file~=t_num)
        disp ('Looks like the calculation is suspended! Please check.')
    end
end