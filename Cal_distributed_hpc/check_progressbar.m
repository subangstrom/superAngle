%check parallel run progress remotely 
clear
txt_file=dir('*.txt');

for i=1:length(txt_file)
    filename=txt_file(i).name;
    A=importdata(filename);
    disp(['File:',filename,' | Job is ',num2str((length(A)-1)/A(1)*100),'% done'])
end
