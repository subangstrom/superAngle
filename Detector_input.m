function [ detector_para, tot_Det_num ] = Detector_input( Detector_file )
% Calculate x-ray with eular angles (in sphere coordinate) that enters
% the EDS detectors
% Detector parameters must be input from file
% Copyright by Weizong Xu, email to wxu4@ncsu.edu for any info and questions. 
% March, 2015

%clear;

%Detector_file = 'detector.xlsx';


% if exist(Detector_file, 'file')
%     Datainput = xlsread(Detector_file);
%     tot_Det_num = max (Datainput(:,1));
%     
%     if (tot_Det_num<1)
%         uiwait(msgbox('No detector input, check file format,, set default for Super-X'));
%         tot_Det_num = 4; %Default setting
%         Datainput = [1,1,13,18,135,0,3.1;2,2,13,18,45,0,3.1;3,3,13,18,315,0,3.1;4,4,13,18,225,0,3.1];
%     end
%     
% else
%     uiwait(msgbox('Unable to find input file for EDS detector, set default for Super-X'));
%     tot_Det_num = 4; %Default setting
%     Datainput = [1,1,13,18,135,0,3.1;2,2,13,18,45,0,3.1;3,3,13,18,315,0,3.1;4,4,13,18,225,0,3.1];
% end




[ Datainput ] = Excel_input( Detector_file );
tot_Det_num = max (Datainput(:,1));
    
    if (tot_Det_num<1)
    uiwait(msgbox('Unable to find input file for EDS detector, set default for Super-X'));
    tot_Det_num = 4; %Default setting
    Datainput = [1,1,12,18,135,0,2.9;2,2,12,18,45,0,2.9;3,3,12,18,315,0,2.9;4,4,12,18,225,0,2.9];
    end




detector_para=Datainput;

end

