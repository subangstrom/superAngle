function [ Datainput ] = Excel_input( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if exist(filename, 'file')
    Datainput = xlsread(filename);
else
    uiwait(msgbox(['Error!! Unable to find input file ', filename, '. Please find it.']));
    [S_filename, S_pathname] = uigetfile({'*.xlsx'},['Find input file: ',filename]);
    if (S_filename ~=0)
        Datainput = xlsread([S_pathname,S_filename]);
        uiwait(msgbox(['Set ', S_pathname, ' as default folder.']));
        cd(S_pathname);
    else
        Datainput =0;
    end
end

end