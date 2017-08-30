function [ output, output_omega ] = single_spot( TiltX, TiltY, Detector, sample_para, holder_para, SpuriousX)
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
    EleA_num = sample_para.EleA_num;
    EleA_shell = sample_para.EleA_shell;
    EleB_num = sample_para.EleB_num;
    EleB_shell = sample_para.EleB_shell;
    
    %Get detector efficiency;
    temp = Excel_input ('Atomic_info.xlsx');
    Xray_table (:,:,1) = temp(1:3,:); % K-shell xray energy, crossection, flurence yield 
    Xray_table (:,:,2) = temp(4:6,:); % L-shell xray energy, crossection, flurence yield 
    Xray_table (:,:,3) = temp(7:9,:); % M-shell xray energy, crossection, flurence yield 
    EleA_energy = Xray_table(1,EleA_num,EleA_shell); %keV
    EleB_energy = Xray_table(1,EleB_num,EleB_shell); %keV
    [ Det_eff_A, Det_eff_B ] = get_Det_eff( EleA_energy, EleB_energy, 'SDD_windowless_efficiency.xlsx' );

    for i=1:tot_Det_num;
        [A_Pout(i), B_Pout(i)] = Point_search(sample_para, holder_para, angle_search(:,:,i), SpuriousX(:,i));
        Pratio(i)= A_Pout(i)/B_Pout(i); %A - Al; B - Ni; 
        output (i,:) = [i, A_Pout(i)/sum(angle_search(:,5,i)),B_Pout(i)/sum(angle_search(:,5,i)),Pratio(i)]; % counts*1e3;
        output_omega (i,:) = [i, A_Pout(i),B_Pout(i),Pratio(i)]; % counts*1e3;
        output_omega_deteff (i,:) = [i, A_Pout(i)*Det_eff_A,B_Pout(i)*Det_eff_B,Pratio(i)*Det_eff_A/Det_eff_B];
    end
    
    
    sym_A = get_element_name(EleA_num,EleA_shell);
    sym_B = get_element_name(EleB_num,EleB_shell);
    output (i+1,:) = [0.1234, sum(A_Pout)/sum(angle_search(:,5,i)),sum(B_Pout)/sum(angle_search(:,5,i)),sum(A_Pout)/sum(B_Pout)];  %all detectors
    output_omega (i+1,:) = [0.1234, sum(A_Pout),sum(B_Pout),sum(A_Pout)/sum(B_Pout)];  %all detectors
    disp(['TiltX: ',num2str(TiltX),' degree   TiltY: ',num2str(TiltY), ' degree']) 
    disp('___________________________________________') 
    disp('Effective Solid angle reaching to detectors')  
    disp(['detector#   ',sym_A,'        ',sym_B, '    ', sym_A,'/', sym_B,' ratio'])
    disp(num2str(output_omega));
    disp('  ')
    disp('___________________________________________')
    output_omega_deteff (i+1,:) = [0.1234, sum(A_Pout)*Det_eff_A,sum(B_Pout)*Det_eff_B,sum(A_Pout)/sum(B_Pout)*Det_eff_A/Det_eff_B];  %all detectors
    disp('Effective Solid angle * detector efficiency')  
    disp(['detector#   ',sym_A,'        ',sym_B, '    ', sym_A,'/', sym_B,' ratio'])
    disp(num2str(output_omega_deteff));

end

