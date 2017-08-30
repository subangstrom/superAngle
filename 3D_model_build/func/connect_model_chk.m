function [ model_connect, chk ] = connect_model_chk( model_sample, model_ether, model_connect )
%Connect specimen model and vacuum model
%check if some tetrahderal is not correctly registered in the neighbour
%list
%Weizong Xu, August, 2017

model_connect_comp=model_connect;
%% must do another check
% must consider this effect for a complete set
%if two-point line from ether/sample intercept with sample/ether tetrahedron,
%neighbour of this line is also with that tetrahedron
l=size(model_connect.tetra_reg_sample,1);
h = waitbar(0, 'Calculating...');
tic
for ii=1:l
    tmpe=model_connect.tetra_reg_sample{ii,1}; %get neighbour tetra from ether
    for jj=1:length(tmpe)
        for kk=1:6
            mesh_num=ii;
            line_num=model_ether.nei_list_coline.line{tmpe(jj),kk};
            line_p=model_ether.p(:,line_num(1));
            line_p=line_p';
            line_direction=model_ether.p(:,line_num(1))-model_ether.p(:,line_num(2));
            line_direction=line_direction';
            [  P1, ~, P2, ~, chk_out ] = intercept_line_tetra_fast( line_p, line_direction, mesh_num, model_sample.p, model_sample.t_sort, model_sample.nei_list_coplane );
            if (chk_out>=1)
                A1=model_ether.p(:,line_num(1))';
                A2=model_ether.p(:,line_num(2))';
                A1A2=sum((A1-A2).^2);
                P1P2=sum((P1-P2).^2);
                A1P1=sum((A1-P1).^2);
                A1P2=sum((A1-P2).^2);
                A2P1=sum((A2-P1).^2);
                A2P2=sum((A2-P2).^2);
                if ((A1P2>A1A2 && A1P2>P1P2 && A1P2>A1P1 && A1P2>A2P1 && A1P2>A2P2) || ...
                        (A2P2>A1A2 && A2P2>P1P2 && A2P2>A1P1 && A2P2>A2P1 && A2P2>A1P2) || ...
                        (A1P1>A1A2 && A1P1>P1P2 && A1P1>A2P2 && A1P1>A2P1 && A1P1>A1P2) || ...
                        (A2P1>A1A2 && A2P1>P1P2 && A2P1>A1P1 && A2P1>A2P2 && A2P1>A1P2))
                    % not intercept with tetrahderon
                else
                    model_connect.tetra_reg_sample{ii,1}=[model_connect.tetra_reg_sample{ii,1} model_ether.nei_list_coline.seperate{tmpe(jj),kk}];
                end
            end
        end
    end
    waitbar(ii/l, h);
end
close (h);
toc

l=size(model_connect.tetra_reg_ether,1);
h = waitbar(0, 'Calculating...');
tic
for ii=1:l
    tmpe=model_connect.tetra_reg_ether{ii,1}; %get neighbour tetra from ether
    for jj=1:length(tmpe)
        for kk=1:6
            mesh_num=ii;
            line_num=model_sample.nei_list_coline.line{tmpe(jj),kk};
            line_p=model_sample.p(:,line_num(1));
            line_p=line_p';
            line_direction=model_sample.p(:,line_num(1))-model_sample.p(:,line_num(2));
            line_direction=line_direction';
            [  P1, ~, P2, ~, chk_out ] = intercept_line_tetra_fast( line_p, line_direction, mesh_num, model_ether.p, model_ether.t_sort, model_ether.nei_list_coplane );
            if (chk_out>=1)
                A1=model_sample.p(:,line_num(1))';
                A2=model_sample.p(:,line_num(2))';
                A1A2=sum((A1-A2).^2);
                P1P2=sum((P1-P2).^2);
                A1P1=sum((A1-P1).^2);
                A1P2=sum((A1-P2).^2);
                A2P1=sum((A2-P1).^2);
                A2P2=sum((A2-P2).^2);
                if ((A1P2>A1A2 && A1P2>P1P2 && A1P2>A1P1 && A1P2>A2P1 && A1P2>A2P2) || ...
                        (A2P2>A1A2 && A2P2>P1P2 && A2P2>A1P1 && A2P2>A2P1 && A2P2>A1P2) || ...
                        (A1P1>A1A2 && A1P1>P1P2 && A1P1>A2P2 && A1P1>A2P1 && A1P1>A1P2) || ...
                        (A2P1>A1A2 && A2P1>P1P2 && A2P1>A1P1 && A2P1>A2P2 && A2P1>A1P2))
                    % not intercept with tetrahderon
                else
                    model_connect.tetra_reg_ether{ii,1}=[model_connect.tetra_reg_ether{ii,1} model_sample.nei_list_coline.seperate{tmpe(jj),kk}];
                end
            end
        end
    end
    waitbar(ii/l, h);
end
close (h);
toc

%% sort data and get rid of duplicated point, AGAIN
for ii=1:size(model_connect.tetra_reg_sample,1)
    model_connect.tetra_reg_sample{ii,1}=unique(model_connect.tetra_reg_sample{ii,1});
end

for ii=1:size(model_connect.tetra_reg_ether,1)
    model_connect.tetra_reg_ether{ii,1}=unique(model_connect.tetra_reg_ether{ii,1});
end
%% sample and ether register with each other, AGAIN
l=size(model_connect.tetra_reg_sample,1);
for ii=1:l
    tmpc=model_connect.tetra_reg_sample{ii,1};
    if (isempty(tmpc)==0)
        for jj=1:length(tmpc)           
            model_connect.tetra_reg_ether{tmpc(jj),1}=[model_connect.tetra_reg_ether{tmpc(jj),1} ii];
        end
    end
end

l=size(model_connect.tetra_reg_ether,1);
for ii=1:l
    tmpd=model_connect.tetra_reg_ether{ii,1};
    if (isempty(tmpd)==0)
        for jj=1:length(tmpd)           
            model_connect.tetra_reg_sample{tmpd(jj),1}=[model_connect.tetra_reg_sample{tmpd(jj),1} ii];
        end
    end
end

%% sort data and get rid of duplicated point, AGAIN, AGAIN
for ii=1:size(model_connect.tetra_reg_sample,1)
    model_connect.tetra_reg_sample{ii,1}=unique(model_connect.tetra_reg_sample{ii,1});
end

for ii=1:size(model_connect.tetra_reg_ether,1)
    model_connect.tetra_reg_ether{ii,1}=unique(model_connect.tetra_reg_ether{ii,1});
end
%% check if input and output is same, i.e. all neighbours are found
%must have this step if size of two model is very different, may take a lot of time
chk=isequal(model_connect,model_connect_comp);

end