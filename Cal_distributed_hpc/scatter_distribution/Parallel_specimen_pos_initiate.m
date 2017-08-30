%% Initiate a broadcast value in dropbox
if exist('control.mat','file')
    prompt = 'Exist control.mat file, replace it with a new one? Y/N [N]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
    end
else
    str = 'Y';
end

if strcmp(str, 'Y')
    X_range=[-2000,2000]/25;X_step=200/25;
    Y_range=[100,4100]/25;Y_step=200/25;
    num=0;
    for sx=X_range(1):X_step:X_range(2)  %21 points
        for sy=Y_range(1):Y_step:Y_range(2) %21 points 21*21*61=26901 points 26901*25s/4*0.8=134505s ~37hrs
            num=num+1;
            control_table(num,1)=sx;
            control_table(num,2)=sy;
            control_table(num,3)=1; % if 1 somebody is working on it
        end
    end
    save('control.mat', 'control_table')
    disp('Job is done. File control.mat is saved.')
else
	disp('Job is cancelled.')
end