function [ output, output_omega, output_counts ] = single_spot_3D_display( TiltX, TiltY, Detector, sample_para, holder_para, SpuriousX, model_sample, model_ether, model_connect, ether_start_p)
%Omega calculation at single spot.
%Weizong Xu, Feb. 2015, wxu4@ncsu.edu

    tot_Det_num=Detector.tot_Det_num;
    angle_search=Detector.angle_search;
    sample_para.TiltX=TiltX;
    sample_para.TiltY=TiltY;
    A_Pout=zeros(tot_Det_num,1);
    B_Pout=zeros(tot_Det_num,1);
    Pratio=zeros(tot_Det_num,1);
    output=zeros(tot_Det_num+1,4);
    output_omega=zeros(tot_Det_num+1,4);
    output_omega_deteff=zeros(tot_Det_num+1,4);
    output_counts=zeros(tot_Det_num+1,4);
    EleA_num = sample_para.EleA_num;
    EleA_shell = sample_para.EleA_shell;
    EleB_num = sample_para.EleB_num;
    EleB_shell = sample_para.EleB_shell;
    if (sample_para.chk_2nd==1)
        disp('WARNING! 2nd interception is set to be enable, but it is disable in current visual version!')
    end
    
    %Get detector efficiency;
    temp = Excel_input ('Atomic_info.xlsx');
    Xray_table (:,:,1) = temp(1:3,:); % K-shell xray energy, crossection, flurence yield 
    Xray_table (:,:,2) = temp(4:6,:); % L-shell xray energy, crossection, flurence yield 
    Xray_table (:,:,3) = temp(7:9,:); % M-shell xray energy, crossection, flurence yield 
    EleA_energy = Xray_table(1,EleA_num,EleA_shell); %keV
    EleB_energy = Xray_table(1,EleB_num,EleB_shell); %keV
    [ Det_eff_A, Det_eff_B ] = get_Det_eff( EleA_energy, EleB_energy, 'SDD_windowless_efficiency.xlsx' );

if (sample_para.t_chk==10)
	control.search_step=500*model_sample.geo_real_ratio; %along y direction in current model setting
	control.pos_acc=1*model_sample.geo_real_ratio;
	control.chk_display=-1;
	[ ether_start_p, ~ ] = thickness_const_model( ether_start_p, model_sample, model_ether, model_connect, sample_para, control );
end

    
    [ ether_start_p, chk_contact ] = tilt_compensate_startpoint( model_sample, model_ether, model_connect, sample_para, ether_start_p, 1 );
if (chk_contact==1) %have interception with sample
    [ model_sample, model_ether ] = tilt_model( model_sample, model_ether, sample_para );
    ether_start_p = sample_normal(ether_start_p,sample_para);
    [ model_slice ] = slice_model( model_sample, model_ether, model_connect, sample_para,ether_start_p, 2, 1);
    %pause(0.2);
    u_value=zeros(length(model_sample.p),4);
    for i=1:tot_Det_num
        sample_para.det_num=i;
        [A_Pout(i), B_Pout(i), u_value] = Point_search_3D_display(sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), model_sample, model_ether, model_connect, model_slice, u_value);
        %[A_Pout(i), B_Pout(i)] = Point_search_3D_fast(sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), model_sample, model_slice);
        %[A_Pout(i), B_Pout(i)] = Point_search_3D_fast_test(sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i), model_sample, model_ether, model_connect, model_slice);
        Pratio(i)= A_Pout(i)/B_Pout(i); %A - Al; B - Ni; 
        output (i,:) = [i, A_Pout(i)/sum(angle_search(:,5,i)),B_Pout(i)/sum(angle_search(:,5,i)),Pratio(i)]; % counts*1e3;
        output_omega (i,:) = [i, A_Pout(i),B_Pout(i),Pratio(i)]; % counts*1e3;
        output_omega_deteff (i,:) = [i, A_Pout(i)*Det_eff_A,B_Pout(i)*Det_eff_B,Pratio(i)*Det_eff_A/Det_eff_B];
    end
    
else
    disp('WARNING! No contact!')
    model_slice.Model_Thickness=0;
    model_slice.Real_Thickness=0;
    u_value=zeros(length(model_sample.p),4);
end
    
v_a = version('-release');
v_b = str2double(v_a(regexp(v_a,'\d')));
if (v_b>=2015)
    
    line1=ether_start_p(1);
    line2=ether_start_p(2);
    line31=ether_start_p(3)-50-model_slice.Model_Thickness/2;
    line32=ether_start_p(3)+50+model_slice.Model_Thickness/2;
    
    figure;pdeplot3D(model_sample.p,model_sample.t);
    hold on;
    plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r')
    hold off;
    
   if (chk_contact==1)
    figure;pdeplot3D(model_sample.p,model_sample.t,'colormapdata',u_value(:,1));
    caxis([-max(u_value(:,1))/3 max(u_value(:,1))]);
    hold on;
    pdeplot3D(model_sample.p,model_sample.t,'FaceAlpha',0);
    title('Ideal beam intensity without specimen absorption');
    hold on;
    plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r');
    hold off;
    %view(-165, 21);
    %print(['model_beam_intensity'],'-dpng','-r800');
    
    
    figure;pdeplot3D(model_sample.p,model_sample.t,'colormapdata',u_value(:,2));
    caxis([-max(u_value(:,2))/3 max(u_value(:,2))]);
    hold on;
    pdeplot3D(model_sample.p,model_sample.t,'FaceAlpha',0);
    title('Element A - beam intensity');
    hold on;
    plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r');
    hold off;
    
    figure;pdeplot3D(model_sample.p,model_sample.t,'colormapdata',u_value(:,3));
    caxis([-max(u_value(:,3))/3 max(u_value(:,3))]);
    hold on;
    pdeplot3D(model_sample.p,model_sample.t,'FaceAlpha',0);
    title('Element B - beam intensity');
    hold on;
    plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r');
    hold off;
    
    for kk=1:size(u_value,1)
        if (u_value(kk,1)>0)
            u_value(kk,2)=u_value(kk,2)./u_value(kk,1);
            u_value(kk,3)=u_value(kk,3)./u_value(kk,1);
        end
    end

    figure;pdeplot3D(model_sample.p,model_sample.t,'colormapdata',u_value(:,2));
    caxis([-max(u_value(:,2))/3 max(u_value(:,2))]);
    title('Element A - Correction factor');
    hold on;
    plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r');
    hold off;
    
    figure;pdeplot3D(model_sample.p,model_sample.t,'colormapdata',u_value(:,3));
    title('Element B - Correction factor');
    caxis([-max(u_value(:,3))/3 max(u_value(:,3))]);
    hold on;
    plot3([line1,line1],[line2,line2],[line31,line32],'LineWidth',2,'color','r');
    hold off;
   end
else
    disp ('WARNING! Need Matlab version 2015a or later to display 3D model!')
end

if (chk_contact==1)
    sym_A = get_element_name(EleA_num,EleA_shell);
    sym_B = get_element_name(EleB_num,EleB_shell);
    output (tot_Det_num+1,:) = [0.1234, sum(A_Pout)/sum(angle_search(:,5,tot_Det_num)),sum(B_Pout)/sum(angle_search(:,5,tot_Det_num)),sum(A_Pout)/sum(B_Pout)];  %all detectors
    output_omega (tot_Det_num+1,:) = [0.1234, sum(A_Pout),sum(B_Pout),sum(A_Pout)/sum(B_Pout)];  %all detectors
    disp(['TiltX: ',num2str(TiltX),' degree   TiltY: ',num2str(TiltY), ' degree']) 
    disp('___________________________________________') 
    disp('Effective Solid angle reaching to detectors')  
    disp(['detector#   ',sym_A,'        ',sym_B, '    ', sym_A,'/', sym_B,' ratio'])
    disp(num2str(output_omega));
    disp('  ')
    disp('___________________________________________')
    output_omega_deteff (tot_Det_num+1,:) = [0.1234, sum(A_Pout)*Det_eff_A,sum(B_Pout)*Det_eff_B,sum(A_Pout)/sum(B_Pout)*Det_eff_A/Det_eff_B];  %all detectors
    disp('Effective Solid angle * detector efficiency')  
    disp(['detector#   ',sym_A,'        ',sym_B, '    ', sym_A,'/', sym_B,' ratio'])
    disp(num2str(output_omega_deteff));
    
%if output counts 
if (sample_para.probe_Ne>1 && sample_para.acquire_time>1)
    [ convert_factor_A,convert_factor_B, ~,~ ] = absolute_scale_factor( sample_para );
    convert_factor_A=convert_factor_A/sample_para.Thickness*model_slice.Real_Thickness;
    convert_factor_B=convert_factor_B/sample_para.Thickness*model_slice.Real_Thickness;
    for i=1:tot_Det_num;
        output_counts (i,:) = [i, A_Pout(i)*convert_factor_A,B_Pout(i)*convert_factor_B,Pratio(i)*convert_factor_A/convert_factor_B];
    end
    output_counts (i+1,:) = [0.1234, sum(A_Pout)*convert_factor_A,sum(B_Pout)*convert_factor_B,sum(A_Pout)/sum(B_Pout)*convert_factor_A/convert_factor_B];  %all detectors
    disp('  ')
    disp(['TiltX: ',num2str(TiltX),' degree   TiltY: ',num2str(TiltY), ' degree']) 
    disp('___________________________________________')
    disp('Counts')  
    disp(['detector#   ',sym_A,'        ',sym_B, '    ', sym_A,'/', sym_B,' ratio'])
    disp(num2str(output_counts));
end

end

end

