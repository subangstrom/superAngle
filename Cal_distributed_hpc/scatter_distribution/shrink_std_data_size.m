clear
load('control.mat')
%load('Tilt_search_Result_std.mat','Tilt_map_std_all', 'Detector', 'search_Deg_2D', 'd_Deg_2D', 'sample_para', 'holder_para', 'SpuriousX')
load('Tilt_search_Result_std_big.mat')
Tilt_map_std_all_in=Tilt_map_std_all;
clear Tilt_map_std_all;

num_count=0;
for i=1:size(control_table,1)
    tmp_x=control_table(i,1);
    tmp_y=control_table(i,2);
    dup_chk=0;
    for j=1:size(Tilt_map_std_all_in.ether_start_p,1)
        coor_temp=Tilt_map_std_all_in.ether_start_p{j,1};
        if (tmp_x==coor_temp(1) && tmp_y==coor_temp(2) && dup_chk==0)
            disp(['Find std data for ',num2str(tmp_x),' ',num2str(tmp_y)])            
            Tilt_map_std_all.Tilt_map_std{i,1}=Tilt_map_std_all_in.Tilt_map_std{j,1};
            Tilt_map_std_all.Thickness{i,1}=Tilt_map_std_all_in.Thickness{j,1};
            Tilt_map_std_all.ether_start_p{i,1}=Tilt_map_std_all_in.ether_start_p{j,1};
            num_count=num_count+1;
            dup_chk=1;
        end
    end
end
Tilt_map_std_all.t_chk=Tilt_map_std_all_in.t_chk;

if (num_count==size(control_table,1))
    disp('Find all points! Data is saved.')
    save('Tilt_search_Result_std.mat','Detector','SpuriousX','Tilt_map_std_all','sample_para','holder_para')
else
    disp('Some data is lost, please check!')
end