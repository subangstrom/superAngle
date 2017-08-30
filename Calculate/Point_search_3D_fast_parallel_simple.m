function [ Al_out, Ni_out ] = Point_search_3D_fast_parallel_simple( sample_para, holder_para, angle_search, Spurious, model_sample, model_slice)
%Calculate X-ray absorption ratio and absorption correction efficient
%Search for one point, only consider one interception
%Weizong Xu, Feb. 2015


chk_2nd=sample_para.chk_2nd; %1 consider 2nd interception (25-100x slow), other not consider
if (chk_2nd==1)
    disp('Error! 2nd interception is enable, but not supported in this funtion. Please have a check!')
    return;
end

u1 = sample_para.uA*1e-7;
u2 = sample_para.uB*1e-7;
n = model_slice.N;

tot_Al=0;
Al_out=0;
tot_Ni=0;
Ni_out=0;
[H_angle_out] = Holder_shadow(sample_para, holder_para, angle_search);
temp_num=size(H_angle_out);   
%u_value=zeros(length(p),4);
p=model_sample.p;%fast data transfer
t_sort=model_sample.t_sort;
four_neighbour=model_sample.nei_list_coplane.four_neighbour;
four_neighbour_mesh_num=zeros(1,4);
tetra_P1=zeros(1,3);
tetra_P2=zeros(1,3);
tetra_P3=zeros(1,3);
tetra_P4=zeros(1,3);
ll=size(model_slice.p,1);
mesh_num_bk=1;
for ii=1:ll
    %ii
    tot_Al_tt=0;
    tot_Ni_tt=0;
    sample_start_p=model_slice.p{ii};
    mesh_num_ii=model_slice.n{ii};
    for jj=1:temp_num(1) %size(H_angle_out,1)
        line_direction=H_angle_out(jj,6:8);
        if (jj>1)
            line_point=sample_start_p;
            tetra_P1(1)=p(1,t_sort(1,mesh_num_bk)); %optimize for faster run
            tetra_P2(1)=p(1,t_sort(2,mesh_num_bk));
            tetra_P3(1)=p(1,t_sort(3,mesh_num_bk));
            tetra_P4(1)=p(1,t_sort(4,mesh_num_bk));
            tetra_P1(2)=p(2,t_sort(1,mesh_num_bk));
            tetra_P2(2)=p(2,t_sort(2,mesh_num_bk));
            tetra_P3(2)=p(2,t_sort(3,mesh_num_bk));
            tetra_P4(2)=p(2,t_sort(4,mesh_num_bk));
            tetra_P1(3)=p(3,t_sort(1,mesh_num_bk));
            tetra_P2(3)=p(3,t_sort(2,mesh_num_bk));
            tetra_P3(3)=p(3,t_sort(3,mesh_num_bk));
            tetra_P4(3)=p(3,t_sort(4,mesh_num_bk));
            four_neighbour_mesh_num(1)=four_neighbour(mesh_num_bk,1);
            four_neighbour_mesh_num(2)=four_neighbour(mesh_num_bk,2);
            four_neighbour_mesh_num(3)=four_neighbour(mesh_num_bk,3);
            four_neighbour_mesh_num(4)=four_neighbour(mesh_num_bk,4);
            [  line_point, n_tmp ] = intercept_line_tetra_fast_forward( line_point, line_direction, tetra_P1, tetra_P2, tetra_P3, tetra_P4, four_neighbour_mesh_num);
 
            if (isnan(n_tmp)==0)
                mesh_num_bk=n_tmp; %must update, bug fix
            	tetra_P1(1)=p(1,t_sort(1,n_tmp)); %optimize for faster run
                tetra_P2(1)=p(1,t_sort(2,n_tmp));
                tetra_P3(1)=p(1,t_sort(3,n_tmp));
                tetra_P4(1)=p(1,t_sort(4,n_tmp));
                tetra_P1(2)=p(2,t_sort(1,n_tmp));
                tetra_P2(2)=p(2,t_sort(2,n_tmp));
                tetra_P3(2)=p(2,t_sort(3,n_tmp));
                tetra_P4(2)=p(2,t_sort(4,n_tmp));
                tetra_P1(3)=p(3,t_sort(1,n_tmp));
                tetra_P2(3)=p(3,t_sort(2,n_tmp));
                tetra_P3(3)=p(3,t_sort(3,n_tmp));
                tetra_P4(3)=p(3,t_sort(4,n_tmp));
                four_neighbour_mesh_num(1)=four_neighbour(n_tmp,1);
                four_neighbour_mesh_num(2)=four_neighbour(n_tmp,2);
                four_neighbour_mesh_num(3)=four_neighbour(n_tmp,3);
                four_neighbour_mesh_num(4)=four_neighbour(n_tmp,4);
                [  line_point, n_tmp ] = intercept_line_tetra_fast_forward( line_point, line_direction, tetra_P1, tetra_P2, tetra_P3, tetra_P4, four_neighbour_mesh_num);
                if (isnan(n_tmp)==0)
                    line_point=[];
                end
            end
            
        end

        if (jj==1 || isempty(line_point))
        %[Xray_dist, ~, ~] = get_intercept_model(sample_start_p,mesh_num_ii,line_direction,model_sample,1); %0 get intercept in both direction; >0 in forward direction <0 backward direction
        %Xray_dist = get_intercept_model_forward(sample_start_p, mesh_num_ii, line_direction, p, t_sort, four_neighbour); %0 get intercept in forward direction only 
            line_point=sample_start_p; %equal to call function:get_intercept_model_forward
            mesh_num=mesh_num_ii;
            while (isnan(mesh_num)==0)
                
                tetra_P1(1)=p(1,t_sort(1,mesh_num)); %optimize for faster run
                tetra_P2(1)=p(1,t_sort(2,mesh_num));
                tetra_P3(1)=p(1,t_sort(3,mesh_num));
                tetra_P4(1)=p(1,t_sort(4,mesh_num));
                tetra_P1(2)=p(2,t_sort(1,mesh_num));
                tetra_P2(2)=p(2,t_sort(2,mesh_num));
                tetra_P3(2)=p(2,t_sort(3,mesh_num));
                tetra_P4(2)=p(2,t_sort(4,mesh_num));
                tetra_P1(3)=p(3,t_sort(1,mesh_num));
                tetra_P2(3)=p(3,t_sort(2,mesh_num));
                tetra_P3(3)=p(3,t_sort(3,mesh_num));
                tetra_P4(3)=p(3,t_sort(4,mesh_num));
                four_neighbour_mesh_num(1)=four_neighbour(mesh_num,1);
                four_neighbour_mesh_num(2)=four_neighbour(mesh_num,2);
                four_neighbour_mesh_num(3)=four_neighbour(mesh_num,3);
                four_neighbour_mesh_num(4)=four_neighbour(mesh_num,4);
                
                [  p_forward, n_forward ] = intercept_line_tetra_fast_forward( line_point, line_direction, tetra_P1, tetra_P2, tetra_P3, tetra_P4, four_neighbour_mesh_num);
                line_point=p_forward;
                mesh_num_bk=mesh_num;
                mesh_num=n_forward;
            end
            
            if (ii==ll && isempty(line_point)) %bug fix, at the bottom side of sample point not intercept sample and exit
                line_point=sample_start_p;
            end
            
            if (isempty(line_point))
                %disp('Error'); %fast search mode is not valid, return to
                %normal slower way
                line_point=sample_start_p;
                mesh_num=mesh_num_ii;
                while (isnan(mesh_num)==0)
                    [  p_forward, n_forward, ~, ~, ~] = intercept_line_tetra_fast( line_point, line_direction, mesh_num, model_sample.p,model_sample.t_sort, model_sample.nei_list_coplane);
                    %[  p_forward, n_forward, ~, ~, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_sample );
                    line_point=p_forward;
                    mesh_num_bk=mesh_num;
                    mesh_num=n_forward;
                end
            end
        end
            if (isempty(line_point))
                disp('Error in calculating distance. Simulation may continue using neighbour ray distance value.')
                disp(['At ii=',num2str(ii),' jj=',num2str(jj)]);
            else
                d_intercept=line_point-sample_start_p;
                Xray_dist=sqrt(d_intercept(1)^2+d_intercept(2)^2+d_intercept(3)^2);
            end
                
        tot_Al_tt=tot_Al_tt+H_angle_out(jj,3)*exp(-(u1*Xray_dist*model_slice.RealModel_ratio));
        tot_Ni_tt=tot_Ni_tt+H_angle_out(jj,4)*exp(-(u2*Xray_dist*model_slice.RealModel_ratio));
%         if (ii==round(size(model_slice.p,1)/2))
%             list_p=model_sample.t_sort(:,mesh_num_bk);
%             for kk=1:4
%                 u_value(list_p(kk),1)=u_value(list_p(kk),1)+1;
%                 u_value(list_p(kk),2)=u_value(list_p(kk),2)+exp(-(u1*Xray_dist*model_slice.RealModel_ratio));%Al
%                 u_value(list_p(kk),3)=u_value(list_p(kk),3)+exp(-(u2*Xray_dist*model_slice.RealModel_ratio));%Ni
%                 u_value(list_p(kk),4)=mesh_num_bk;
%             end
%         end
    end
    tot_Al=tot_Al+tot_Al_tt;
    tot_Ni=tot_Ni+tot_Ni_tt;
end

%figure;pdeplot3D(p,model_sample.t,'colormapdata',u_value(:,1));

if (temp_num(1)>=1)
    Al_out=tot_Al/n+Spurious(1);  %avg Absorption Al
    Ni_out=tot_Ni/n+Spurious(2);  %avg Absorption Ni
else
    Al_out=0;
    Ni_out=1e-20;
end

end