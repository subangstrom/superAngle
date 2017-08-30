function [ Al_out, Ni_out ] = Point_search_3D_fast_test( sample_para, holder_para, angle_search, Spurious, model_sample, model_ether, model_connect, model_slice)
%Calculate X-ray absorption ratio and absorption correction efficient
%Search for one point
%Weizong Xu, Feb. 2015

chk_2nd=sample_para.chk_2nd; %1 consider 2nd interception (25-100x slow), other not consider

u1 = sample_para.uA*1e-7;
u2 = sample_para.uB*1e-7;
n = model_slice.N;

tot_Al=0;
Al_out=0;
tot_Ni=0;
Ni_out=0;

[H_angle_out] = Holder_shadow(sample_para, holder_para, angle_search);
temp_num=size(H_angle_out);   
p=model_sample.p;%fast data transfer
t_sort=model_sample.t_sort;
e_tetra_volume=model_ether.tetra_volume;
tetra_volume=model_sample.tetra_volume;
e_p=model_ether.p;
e_t_sort=model_ether.t_sort;
four_neighbour=model_sample.nei_list_coplane.four_neighbour;
e_four_neighbour=model_ether.nei_list_coplane.four_neighbour;
four_neighbour_mesh_num=zeros(1,4);
tetra_P1=zeros(1,3);
tetra_P2=zeros(1,3);
tetra_P3=zeros(1,3);
tetra_P4=zeros(1,3);
ll=size(model_slice.p,1);
% mesh_num_bk_2nd=-1;
mesh_num_bk=1;
%h = waitbar(0, ['Working on detector #', num2str(sample_para.det_num)]);
%ppp=ProgressBar(ll);
for ii=1:ll
    %ii
    tot_Al_tt=0;
    tot_Ni_tt=0;
    sample_start_p=model_slice.p{ii};
    mesh_num_ii=model_slice.n{ii};
    for jj=1:temp_num(1) %size(H_angle_out,1)
        line_direction=H_angle_out(jj,6:8);
        if (jj>1) %pretest using last point result, in most cases this will speed up
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
 
            if (isnan(n_tmp)==0) %inside the sample, search it again
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
        
      if (jj==1 || isempty(line_point)) %if not, do regular search
            line_point=sample_start_p; %equal to call function:get_intercept_model_forward
            mesh_num=mesh_num_ii;
            while (isnan(mesh_num)==0)
                ta=t_sort(1,mesh_num); tb=t_sort(2,mesh_num); tc=t_sort(3,mesh_num); td=t_sort(4,mesh_num);                
                tetra_P1(1)=p(1,ta); tetra_P1(2)=p(2,ta); tetra_P1(3)=p(3,ta);
                tetra_P2(1)=p(1,tb); tetra_P2(2)=p(2,tb); tetra_P2(3)=p(3,tb);
                tetra_P3(1)=p(1,tc); tetra_P3(2)=p(2,tc); tetra_P3(3)=p(3,tc);
                tetra_P4(1)=p(1,td); tetra_P4(2)=p(2,td); tetra_P4(3)=p(3,td);
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
      end
            %***deal with double intercept situation***
            if (chk_2nd==1 && ~isempty(line_point)) %add check, bug fix if line_point is empty.
%             chk_2nd_search=0;
%             %now check the ether point of line_point 
             e_line_point=line_point;
%             if ( mesh_num_bk_2nd==mesh_num_bk )
%                 ether_num=reg_point_tetrahedron_fast2( ether_list_bk_2nd, e_line_point, e_p, e_t_sort, e_tetra_volume );
%                 if (ether_num>0)
%                     chk_2nd_search=1;
%                 end
%             end
                
%             if (chk_2nd_search==0)
            ether_list=model_connect.tetra_reg_sample{mesh_num_bk};
            if (length(ether_list)>1)
                ether_num=reg_point_tetrahedron_fast2( ether_list, e_line_point, e_p, e_t_sort, e_tetra_volume );
%                 ether_list_bk_2nd=ether_num;
%                 mesh_num_bk_2nd=mesh_num_bk;
            else
                ether_num=ether_list;
            end
%             end
            
            %find its neighbour tetrahedral along line_direction, i.e. skip
            %this tetrahedral to speed up simulation
            sample_num=0;  
            distance_2nd=0;
            while (sample_num==0 && isnan(ether_num)==0 && distance_2nd==0 && ether_num~=0)
            %get intercept point of this ether (next ether)
                ta=e_t_sort(1,ether_num); tb=e_t_sort(2,ether_num); tc=e_t_sort(3,ether_num); td=e_t_sort(4,ether_num);                
                tetra_P1(1)=e_p(1,ta); tetra_P1(2)=e_p(2,ta); tetra_P1(3)=e_p(3,ta);
                tetra_P2(1)=e_p(1,tb); tetra_P2(2)=e_p(2,tb); tetra_P2(3)=e_p(3,tb);
                tetra_P3(1)=e_p(1,tc); tetra_P3(2)=e_p(2,tc); tetra_P3(3)=e_p(3,tc);
                tetra_P4(1)=e_p(1,td); tetra_P4(2)=e_p(2,td); tetra_P4(3)=e_p(3,td);
                four_neighbour_mesh_num(1)=e_four_neighbour(ether_num,1);
                four_neighbour_mesh_num(2)=e_four_neighbour(ether_num,2);
                four_neighbour_mesh_num(3)=e_four_neighbour(ether_num,3);
                four_neighbour_mesh_num(4)=e_four_neighbour(ether_num,4);               
                [  p_forward, n_forward ] = intercept_line_tetra_fast_forward( e_line_point, line_direction, tetra_P1, tetra_P2, tetra_P3, tetra_P4, four_neighbour_mesh_num);
                if isempty(p_forward) %very rare case
                    disp('Outside tetrahedral. Something is wrong!')
                end
                if (isnan(n_forward)==0)
                    %check if any sample tetrahedron within this ether
                    sample_list=model_connect.tetra_reg_ether{n_forward};
                    s_line_point=p_forward;
                    if (length(sample_list)>=1)
                        sample_num=reg_point_tetrahedron_fast2( sample_list, s_line_point, p, t_sort, tetra_volume );
                        if (sample_num>0)
                            %disp('2nd intercept')
                                line_point_2nd=s_line_point;
                                mesh_num_2nd=sample_num;
                                while (isnan(mesh_num_2nd)==0)
                                	ta=t_sort(1,mesh_num_2nd); tb=t_sort(2,mesh_num_2nd); tc=t_sort(3,mesh_num_2nd); td=t_sort(4,mesh_num_2nd);                
                                    tetra_P1(1)=p(1,ta); tetra_P1(2)=p(2,ta); tetra_P1(3)=p(3,ta);
                                    tetra_P2(1)=p(1,tb); tetra_P2(2)=p(2,tb); tetra_P2(3)=p(3,tb);
                                    tetra_P3(1)=p(1,tc); tetra_P3(2)=p(2,tc); tetra_P3(3)=p(3,tc);
                                    tetra_P4(1)=p(1,td); tetra_P4(2)=p(2,td); tetra_P4(3)=p(3,td);
                                    four_neighbour_mesh_num(1)=four_neighbour(mesh_num_2nd,1);
                                    four_neighbour_mesh_num(2)=four_neighbour(mesh_num_2nd,2);
                                    four_neighbour_mesh_num(3)=four_neighbour(mesh_num_2nd,3);
                                    four_neighbour_mesh_num(4)=four_neighbour(mesh_num_2nd,4);               
                                    [  p_forward, n_forward ] = intercept_line_tetra_fast_forward( line_point_2nd, line_direction, tetra_P1, tetra_P2, tetra_P3, tetra_P4, four_neighbour_mesh_num);
                                    %[  p_forward, n_forward, ~, ~, ~ ] = intercept_line_tetra( line_point_2nd, line_direction, mesh_num_2nd, model_sample );
                                    line_point_2nd=p_forward;
                                    mesh_num_2nd=n_forward;
                                end
                                p_end=line_point_2nd;

                                line_point_2nd=s_line_point;
                                mesh_num_2nd=sample_num;
                                while (isnan(mesh_num)==0)
                                	ta=t_sort(1,mesh_num_2nd); tb=t_sort(2,mesh_num_2nd); tc=t_sort(3,mesh_num_2nd); td=t_sort(4,mesh_num_2nd);                
                                    tetra_P1(1)=p(1,ta); tetra_P1(2)=p(2,ta); tetra_P1(3)=p(3,ta);
                                    tetra_P2(1)=p(1,tb); tetra_P2(2)=p(2,tb); tetra_P2(3)=p(3,tb);
                                    tetra_P3(1)=p(1,tc); tetra_P3(2)=p(2,tc); tetra_P3(3)=p(3,tc);
                                    tetra_P4(1)=p(1,td); tetra_P4(2)=p(2,td); tetra_P4(3)=p(3,td);
                                    four_neighbour_mesh_num(1)=four_neighbour(mesh_num_2nd,1);
                                    four_neighbour_mesh_num(2)=four_neighbour(mesh_num_2nd,2);
                                    four_neighbour_mesh_num(3)=four_neighbour(mesh_num_2nd,3);
                                    four_neighbour_mesh_num(4)=four_neighbour(mesh_num_2nd,4);               
                                    [  p_backward, n_backward ] = intercept_line_tetra_fast_forward( line_point_2nd, -line_direction, tetra_P1, tetra_P2, tetra_P3, tetra_P4, four_neighbour_mesh_num);                                    
                                    %[  ~, ~, p_backward, n_backward, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_sample );
                                    line_point_2nd=p_backward;
                                    mesh_num_2nd=n_backward;
                                end
                                    %p_start=line_point_2nd;
                                    p_intercept=p_end-line_point_2nd;
                                    distance_2nd=sqrt(p_intercept(1)^2+p_intercept(2)^2+p_intercept(3)^2);
                        end
                    else
                        sample_num=0;
                    end
                end
                e_line_point=p_forward;
                ether_num=n_forward;
            end
%            disp('OK')
           else
               distance_2nd=0;
           end        
            %disp(num2str(ether_num));
% %             if (isempty(line_point)) %bug fix if fast intercept_line_tetra failed
% %                 %the bug may happen if point is out of tetrahedra, i.e. the accuracy is not high enough
% %                 line_point=sample_start_p;
% %                 mesh_num=mesh_num_ii;
% %                 while (isnan(mesh_num)==0)
% %                     [  p_forward, n_forward, ~, ~, ~] = intercept_line_tetra_fast( line_point, line_direction, mesh_num, model_sample.p,model_sample.t_sort, model_sample.nei_list_coplane);
% %                     %[  p_forward, n_forward, ~, ~, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_sample );
% %                     line_point=p_forward;
% %                     mesh_num_bk=mesh_num;
% %                     mesh_num=n_forward;
% %                 end
% %             end
%       end
            if (isempty(line_point)) %if bug issue is not solved, try to use its neighbour value.
                tic
                disp('Error in calculating distance. Simulation may continue using neighbour ray distance value.')
                disp(['At ii=',num2str(ii),' jj=',num2str(jj)]);
                toc
            else
                d_intercept=line_point-sample_start_p;
                Xray_dist=sqrt(d_intercept(1)^2+d_intercept(2)^2+d_intercept(3)^2);
           end
                
        tot_Al_tt=tot_Al_tt+H_angle_out(jj,3)*exp(-(u1*(Xray_dist+distance_2nd)*model_slice.RealModel_ratio));
        tot_Ni_tt=tot_Ni_tt+H_angle_out(jj,4)*exp(-(u2*(Xray_dist+distance_2nd)*model_slice.RealModel_ratio));
    end
    tot_Al=tot_Al+tot_Al_tt;
    tot_Ni=tot_Ni+tot_Ni_tt;
    %waitbar(ii/ll, h);
    %ppp.progress;
end
%ppp.stop;
%figure;pdeplot3D(p,model_sample.t,'colormapdata',u_value(:,1));

if (temp_num(1)>=1)
    Al_out=tot_Al/n+Spurious(1);  %avg Absorption Al
    Ni_out=tot_Ni/n+Spurious(2);  %avg Absorption Ni
else
    Al_out=0;
    Ni_out=1e-20;
end
%close(h)
end