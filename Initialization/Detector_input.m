function [ Detector ] = Detector_input( Detector_file, dAngle )
% Calculate x-ray with eular angles (in sphere coordinate) that enters
% the EDS detectors
% Detector parameters must be input from file
% Copyright by Weizong Xu, email to wxu4@ncsu.edu for any info and questions. 
% June, 2016




[ Datainput ] = Excel_input( Detector_file );
tot_Det_num = max (Datainput(:,1));
    
if (tot_Det_num<1)
    uiwait(msgbox('Unable to find input file for EDS detector, set default for Super-X'));
    tot_Det_num = 4; %Default setting
    Datainput = [1,1,12,18,135,0,2.9;2,2,12,18,45,0,2.9;3,3,12,18,315,0,2.9;4,4,12,18,225,0,2.9];
end

% Determine if we're using tables or datasets:
isBeforeR2013b = verLessThan('matlab', '8.2');
if isBeforeR2013b
	detector_para = dataset(Datainput(:,1), Datainput(:,3), Datainput(:,4), Datainput(:,5), Datainput(:,6), Datainput(:,7), Datainput(:,1)*nan, 'VarNames', ...
                {'Detector', 'Distance', 'Takeoff_angle', 'Azimuth_angle', 'Self_tilt', 'Detector_radius', 'Solid_angle'});
else
	detector_para = table(Datainput(:,1), Datainput(:,3), Datainput(:,4), Datainput(:,5), Datainput(:,6), Datainput(:,7), Datainput(:,1)*nan, 'VariableNames', ...
                {'Detector', 'Distance', 'Takeoff_angle', 'Azimuth_angle', 'Self_tilt', 'Detector_radius', 'Solid_angle'});
end

Detector.detector_para=detector_para;
Detector.tot_Det_num=tot_Det_num;
Detector.dAngle=dAngle;
[ Detector] = Detector_setup( Detector );
%detector_para=Datainput;

end

