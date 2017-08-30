function [ A_Pout, B_Pout, A_Pout_counts, B_Pout_counts, p_tetra_num ] = single_spot_3D_parallel( TiltX, TiltY, Detector, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect,ether_start_p,p_tetra_num)
%Omega calculation at single spot.
%Weizong Xu, Feb. 2015, wxu4@ncsu.edu
% tic
    chk_2nd=sample_para.chk_2nd; %1 consider 2nd interception (25-100x slow), other not consider
    tot_Det_num=Detector.tot_Det_num;
    angle_search=Detector.angle_search;
    sample_para.TiltX=TiltX;
    sample_para.TiltY=TiltY;
    sample_para.p_tetra_num=p_tetra_num;
    A_Pout=zeros(tot_Det_num,1);
    B_Pout=zeros(tot_Det_num,1);
    A_Pout_counts=zeros(tot_Det_num,1);
    B_Pout_counts=zeros(tot_Det_num,1);
%     Pratio=zeros(tot_Det_num,1);
%     output=zeros(tot_Det_num+1,4);
%     output_omega_A=zeros(tot_Det_num,2);
%     output_omega_deteff=zeros(tot_Det_num+1,4);
%     EleA_num = sample_para.EleA_num;
%     EleA_shell = sample_para.EleA_shell;
%     EleB_num = sample_para.EleB_num;
%     EleB_shell = sample_para.EleB_shell;
    
%     %Get detector efficiency;
%     temp = Excel_input ('Atomic_info.xlsx');
%     Xray_table (:,:,1) = temp(1:3,:); % K-shell xray energy, crossection, flurence yield 
%     Xray_table (:,:,2) = temp(4:6,:); % L-shell xray energy, crossection, flurence yield 
%     Xray_table (:,:,3) = temp(7:9,:); % M-shell xray energy, crossection, flurence yield 
%     EleA_energy = Xray_table(1,EleA_num,EleA_shell); %keV
%     EleB_energy = Xray_table(1,EleB_num,EleB_shell); %keV
%     [ Det_eff_A, Det_eff_B ] = get_Det_eff( EleA_energy, EleB_energy, 'SDD_windowless_efficiency.xlsx' );
  
[ ether_start_p, chk_contact ] = tilt_compensate_startpoint( model_sample, model_ether, model_connect, sample_para, ether_start_p, 1 );

if (chk_contact==1) %have interception with sample
    [ model_sample, model_ether ] = tilt_model( model_sample, model_ether, sample_para );
    ether_start_p = sample_normal(ether_start_p,sample_para);
    [ model_slice, sample_para ] = slice_model( model_sample, model_ether, model_connect, sample_para,ether_start_p, 0, 1);
    p_tetra_num=sample_para.p_tetra_num;
    %pause(0.2);
    %u_value=zeros(length(model_sample.p),4);
    for i=1:tot_Det_num;
        %[A_Pout(i), B_Pout(i), u_value] = Point_search_3D(sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), model_sample, model_ether, model_connect, model_slice, u_value);
        %[A_Pout(i), B_Pout(i)] = Point_search_3D_fast(sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), model_sample, model_slice);
        sample_para.det_num=i;
        if (chk_2nd==1)
            [A_Pout(i), B_Pout(i)] = Point_search_3D_fast_test(sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), model_sample, model_ether, model_connect, model_slice);
        else
            [A_Pout(i), B_Pout(i)] = Point_search_3D_fast_parallel_simple( sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), model_sample, model_slice);
        end
%         if (isequal(A_Pout(i),A_Pout1(i))==0 || isequal(B_Pout(i),B_Pout1(i))==0)
%            disp('Not consistent') 
%         else
%             disp('equal')
%         end

        %Pratio(i)= A_Pout(i)/B_Pout(i); %A - Al; B - Ni; 
        %output (i,:) = [i, A_Pout(i)/sum(angle_search(:,5,i)),B_Pout(i)/sum(angle_search(:,5,i)),Pratio(i)]; % counts*1e3;
%        output_omega (i,:) = [A_Pout(i),B_Pout(i)]; % counts*1e3;
        %output_omega_deteff (i,:) = [i, A_Pout(i)*Det_eff_A,B_Pout(i)*Det_eff_B,Pratio(i)*Det_eff_A/Det_eff_B];
    end
    
%     line1=ether_start_p(1);
%     line2=ether_start_p(2);
%     line31=ether_start_p(3)-20;
%     line32=ether_start_p(3)+20;
%     
%     figure;pdeplot3D(model_sample.p,model_sample.t);
%     hold on;
%     plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r')
%     hold off;
%     for kk=1:size(u_value,1)
%         if (u_value(kk,1)>0)
%             u_value(kk,2)=u_value(kk,2)./u_value(kk,1);
%             u_value(kk,3)=u_value(kk,3)./u_value(kk,1);
%         end
%     end
%     figure;pdeplot3D(model_sample.p,model_sample.t,'colormapdata',u_value(:,1));
%     caxis([-max(u_value(:,1))/3 max(u_value(:,1))])
%     hold on;
%     pdeplot3D(model_sample.p,model_sample.t,'FaceAlpha',0)
%     hold on;
%     plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r')
%     hold off;
%     
%     figure;pdeplot3D(model_sample.p,model_sample.t,'colormapdata',u_value(:,2));
%     caxis([-max(u_value(:,2))/3 max(u_value(:,2))])
%     hold on;
%     plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r')
%     hold off;
%     
%     figure;pdeplot3D(model_sample.p,model_sample.t,'colormapdata',u_value(:,3));
%     caxis([-max(u_value(:,3))/3 max(u_value(:,3))])
%     hold on;
%     plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r')
%     hold off;  
    
%     sym_A = get_element_name(EleA_num,EleA_shell);
%     sym_B = get_element_name(EleB_num,EleB_shell);
    %output (tot_Det_num+1,:) = [0.1234, sum(A_Pout)/sum(angle_search(:,5,tot_Det_num)),sum(B_Pout)/sum(angle_search(:,5,tot_Det_num)),sum(A_Pout)/sum(B_Pout)];  %all detectors
%     output_omega (tot_Det_num+1,:) = [0.1234, sum(A_Pout),sum(B_Pout),sum(A_Pout)/sum(B_Pout)];  %all detectors
%     disp(['TiltX: ',num2str(TiltX),' degree   TiltY: ',num2str(TiltY), ' degree']) 
%     disp('___________________________________________') 
%     disp('Effective Solid angle reaching to detectors')  
%     disp(['detector#   ',sym_A,'        ',sym_B, '    ', sym_A,'/', sym_B,' ratio'])
%     disp(num2str(output_omega));
%     disp('  ')
%     disp('___________________________________________')
%     output_omega_deteff (tot_Det_num+1,:) = [0.1234, sum(A_Pout)*Det_eff_A,sum(B_Pout)*Det_eff_B,sum(A_Pout)/sum(B_Pout)*Det_eff_A/Det_eff_B];  %all detectors
%     disp('Effective Solid angle * detector efficiency')  
%     disp(['detector#   ',sym_A,'        ',sym_B, '    ', sym_A,'/', sym_B,' ratio'])
%     disp(num2str(output_omega_deteff));
% toc

%if output counts 
if (sample_para.probe_Ne>1 && sample_para.acquire_time>1)
    %[ convert_factor_A,convert_factor_B, ~,~ ] = absolute_scale_factor( sample_para );
    convert_factor_A=sample_para.convert_factor_A;
    convert_factor_B=sample_para.convert_factor_B;
    if (sample_para.t_chk==0)
        convert_factor_A=convert_factor_A/sample_para.Thickness*model_slice.Model_Thickness;
        convert_factor_B=convert_factor_B/sample_para.Thickness*model_slice.Model_Thickness;
    end
    A_Pout_counts=A_Pout*convert_factor_A;
    B_Pout_counts=B_Pout*convert_factor_B;
    %     for i=1:tot_Det_num;
%         output_counts (i,:) = [i, A_Pout(i)*convert_factor_A,B_Pout(i)*convert_factor_B,Pratio(i)*convert_factor_A/convert_factor_B];
%     end
    
end


end


end

