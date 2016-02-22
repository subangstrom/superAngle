function [ A_exp, A_exp_norm, max_A_exp, B_exp, B_exp_norm, max_B_exp, ratio_exp, ratio_exp_all, tot_Det_num_in] = read_exp_data( exp_file )
%Input experiment counts data
%Weizong Xu, March, 2015, wxu4@ncsu.edu

%Data input
if exist(exp_file, 'file')
    Datainput = xlsread(exp_file);
    temp=size(Datainput);
    tot_Det_num_in=(temp(2)-1)/2;
    
    for i=1:tot_Det_num_in

        A_exp_norm(:,1,i) = Datainput(:,1);
        A_exp_norm(:,2,i) = Datainput(:,1+i);
        A_exp(:,1,i) = Datainput(:,1);
        A_exp(:,2,i) = Datainput(:,1+i);

        B_exp_norm(:,1,i) = Datainput(:,1);
        B_exp_norm(:,2,i) = Datainput(:,1+tot_Det_num_in+i);
        B_exp(:,1,i) = Datainput(:,1);
        B_exp(:,2,i) = Datainput(:,1+tot_Det_num_in+i);        
    end
      
        max_A_exp=max(max(A_exp_norm(:,2,:)));
        max_B_exp=mean(max(B_exp_norm(:,2,:)));
    for i=1:tot_Det_num_in        
        A_exp_norm(:,2,i) = A_exp_norm(:,2,i)/max_A_exp;
        B_exp_norm(:,2,i) = B_exp_norm(:,2,i)/max_B_exp;
    end
    
    
    for i=1:tot_Det_num_in
        ratio_exp(:,1,i) = Datainput(:,1);
        ratio_exp(:,2,i) = Datainput(:,1+i)./Datainput(:,1+tot_Det_num_in+i);
    end
    

    ratio_exp_all(:,1)=Datainput(:,1);
    A_exp_all=A_exp(:,2,1);
    B_exp_all=B_exp(:,2,1);
    for i=2:tot_Det_num_in    
        A_exp_all=A_exp_all+A_exp(:,2,i);
        B_exp_all=B_exp_all+B_exp(:,2,i);
    end
    ratio_exp_all(:,2)=A_exp_all./B_exp_all;
    
    chk = 0;
    
    
    
else
    msgbox('No experiment data file input');
    fprintf('No experiment data file input.\n');
    A_exp = 0;
    A_exp_norm = 0;
    max_A_exp = 0;
    B_exp = 0;
    B_exp_norm = 0;
    max_B_exp = 0;
    ratio_exp = 0;
    ratio_exp_all = 0;
    tot_Det_num_in = 0;
end

end

