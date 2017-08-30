function [ model_connect ] = connect_model( model_sample, model_ether )
%Connect specimen model and vacuum model
%Weizong Xu, August, 2017

%% first check if a point from sample is in ether
l=size(model_sample.p,2);
parfor_progress(l);
tmpa=zeros(l,1);
tmpp=model_sample.p;
tic
%h = waitbar(0, 'Calculating...');
parfor ii=1:l
    point=tmpp(:,ii)';
    [~,p_tetra_num] = reg_point_tetrahedron( point, model_ether, -1 );
    tmpa(ii)=p_tetra_num;
    parfor_progress;
    %waitbar(ii/l, h);
end
parfor_progress(0);
%close (h);
toc

tmp_an=zeros(size(model_ether.t,2),1);
for ii=1:l
    if (isnan(tmpa(ii))==0)
        tmp_an(tmpa(ii))=tmp_an(tmpa(ii))+1;
        if (tmp_an(tmpa(ii))==1)
            model_connect.tetra_reg_ether{tmpa(ii),1}=model_sample.p_nei_tetra{ii,1};%check if tetra%p_tetra_num contains tetra from sample, and their num#
        else
            model_connect.tetra_reg_ether{tmpa(ii),1}=[model_connect.tetra_reg_ether{tmpa(ii),1} model_sample.p_nei_tetra{ii,1}];
        end
    end
end

%% then theck if a point from ether is in sample
l=size(model_ether.p,2);
parfor_progress(l);
tmpb=zeros(l,1);
tmpp=model_ether.p;
tic
%h = waitbar(0, 'Calculating...');
parfor ii=1:l
    point=tmpp(:,ii)';
    [~,p_tetra_num] = reg_point_tetrahedron( point, model_sample, -1 );
    tmpb(ii)=p_tetra_num;
    parfor_progress;
    %waitbar(ii/l, h);
end
parfor_progress(0);
%close (h);
toc

tmp_bn=zeros(size(model_sample.t,2),1);
for ii=1:l
    if (isnan(tmpb(ii))==0)
        tmp_bn(tmpb(ii))=tmp_bn(tmpb(ii))+1;
        if (tmp_bn(tmpb(ii))==1)
            model_connect.tetra_reg_sample{tmpb(ii),1}=model_ether.p_nei_tetra{ii,1};%check if tetra%p_tetra_num contains tetra from ether, and their num#
        else
            model_connect.tetra_reg_sample{tmpb(ii),1}=[model_connect.tetra_reg_sample{tmpb(ii),1} model_ether.p_nei_tetra{ii,1}];%check if tetra%p_tetra_num contains tetra from ether, and their num#
        end
    end
end

%% sample and ether register with each other
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

%% sort data and get rid of duplicated point
for ii=1:size(model_connect.tetra_reg_sample,1)
    model_connect.tetra_reg_sample{ii,1}=unique(model_connect.tetra_reg_sample{ii,1});
end

for ii=1:size(model_connect.tetra_reg_ether,1)
    model_connect.tetra_reg_ether{ii,1}=unique(model_connect.tetra_reg_ether{ii,1});
end

%% must do another check
% must consider this effect for a complete set
%if two-point line from ether/sample intercept with sample/ether tetrahedron,
%neighbour of this line is also with that tetrahedron
%loop until all sample is found
chk=0;iter_num=0;
while (chk==0)
    iter_num=iter_num+1;
    disp(['Iteration number: ', num2str(iter_num)])
    [ model_connect, chk ] = connect_model_chk( model_sample, model_ether, model_connect );
end

%% finish
disp('Done. Two models are now connected.')
end

