function [ model_connect ] = connect_model2( model_sample, model_ether, model_connect )
%Connect specimen model and vacuum model
%step 2, register two sample model
%Weizong Xu, August, 2017

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

