function [ cwText ] = GetCommandWindow()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
disp('  ')
jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
jCmdWin = jDesktop.getClient('Command Window');
jTextArea = jCmdWin.getComponent(0).getViewport.getView;
cwText = char(jTextArea.getText);

%s = num2str((1:1000)');
% h = figure(1);
% uicontrol('Parent',h,...
%           'Units','normalized',...
%           'Position',[0.1,0.1,0.8,0.8],...
%           'HorizontalAlignment','left',...
%           'Style','edit',...
%           'Max',100,...
%           'Enable','inactive',...
%           'String',cwText)

% if (length(cwText)>500)
%     cwText_out = cwText((length(cwText)-500):length(cwText));
% else
%     cwText_out=cwText;
% end
end

