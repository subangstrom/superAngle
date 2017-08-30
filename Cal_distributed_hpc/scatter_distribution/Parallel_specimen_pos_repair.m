%% Repair control.mat
clear;
%% function []=repair_check_control(Data_cell, single spot)
prompt = 'Repair All data or just Single one or fix Conflict? A/S/C [A]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'C';
end

if strcmp(str, 'A')
    load('control.mat', 'control_table')
    for i=1:size(control_table,1)
        control_table(i,3)=1;
    end
    
    for i=1:length(control_table)
        sx=control_table(i,1);
        sy=control_table(i,2);
        filename=['Tilt_search_Result_',num2str(sx),'_',num2str(sy),'.mat'];

        if exist(filename, 'file')
            disp([num2str(sx),'_',num2str(sy),' has been calculated.'])
            control_table(i,3)=0;
        end
    end
    save ('control.mat','control_table');
    disp('control.mat has been reset! New file is saved')
end

if strcmp(str, 'S')
    prompt = 'What is the sx value? ';
    ix = input(prompt);
    disp(['sx = ',num2str(ix)]);
    prompt = 'What is the sy value? ';
    iy = input(prompt);
    disp(['sy = ',num2str(iy)]);
    load('control.mat', 'control_table')
    chk_upd=0;
    for i=1:size(control_table,1)
        sx=control_table(i,1);
        sy=control_table(i,2);
        if (sx==ix && sy==iy)
            if (control_table(i,3)==2);
                control_table(i,3)=1;
                chk_upd=1;
                disp(['Point ',num2str(sx),' ',num2str(sy),' has been reset to not calculated.'])
            else
                disp('Looks like this point does not need to be updated. (chk~=2) Nothing is done')
                chk_upd=-1;
            end
        end
    end
    if (chk_upd==0)
        disp('Could not find this point from database! Nothing is done.')
    else
        if (chk_upd==1)
            save ('control.mat','control_table');
            disp('control.mat has been reset! New file is saved')
        end
    end
end

if strcmp(str, 'C')
    load('control.mat', 'control_table')
    
    for i=1:length(control_table)
        sx=control_table(i,1);
        sy=control_table(i,2);
        filename=['Tilt_search_Result_',num2str(sx),'_',num2str(sy),'.mat'];

        if exist(filename, 'file')
            disp([num2str(sx),'_',num2str(sy),' has been calculated.'])
            control_table(i,3)=0;
        end
    end
    save ('control.mat','control_table');
    disp('control.mat has been reset! New file is saved')
end

if ((strcmp(str, 'A')==0) && (strcmp(str, 'S')==0) && (strcmp(str, 'C')==0))
    disp('Incorrect input! Nothing is done.')
end