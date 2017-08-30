function [model_mesh_nei]=cal_neighbour_list(model_mesh, chk)
%generate neighbour list (4 faces) of each tetrahedral 
%Weizong Xu, August, 2017

if (chk==1 || chk==2)
    
model_mesh_nei=model_mesh;
model_mesh_nei.t_sort=model_mesh.t(1:4,:);
for i=1:size(model_mesh_nei.t_sort,2)
    model_mesh_nei.t_sort(:,i)=sort(model_mesh_nei.t_sort(:,i)); %sort is needed before use ismembc
end

tic
l=size(model_mesh_nei.t_sort,2);
%l=1000;
h = waitbar(0, 'Calculating...');
%model_mesh_nei.nei_list_copoint.all=cell(1,1);
for ii=1:l
    t_num=0;t_num2=0;t_num3=0;
    for jj=1:size(model_mesh_nei.t_sort,2)
%         temp1=intersect(model_mesh_nei.t_sort(:,ii), model_mesh_nei.t_sort(:,jj)); %too slow
%         temp=length(temp1);
        temp1=ismembc(model_mesh_nei.t_sort(:,ii), model_mesh_nei.t_sort(:,jj));
        temp=sum(temp1);
        %note: [0 1 1 1 1 1]=ismembc([1 2 3 4 5 17], [2 3 4 5 7 9 11 17])
        if (temp==3) %find one co-plane neighbour
            t_num=t_num+1;
            if (t_num==1)
                model_mesh_nei.nei_list_coplane.all{ii,1}=jj;
            else
                model_mesh_nei.nei_list_coplane.all{ii,1}=[model_mesh_nei.nei_list_coplane.all{ii,1} jj];
            end
            nei_pts_coplane{ii,t_num}=temp1;
        end
        
        if (temp==2) %find one co-line neighbour, a co-plane neighbour is also a co-line neighbour
            t_num2=t_num2+1;
            if (t_num2==1)
            model_mesh_nei.nei_list_coline.all{ii,1}=jj;
            else
                model_mesh_nei.nei_list_coline.all{ii,1}=[model_mesh_nei.nei_list_coline.all{ii,1} jj];
            end
            nei_pts_coline_t{ii,t_num2}=temp1;
        end
        
        if (temp==1) %find one co-point neighbour, a co-point neighbour is also a co-line or co-plane neighbour
            t_num3=t_num3+1;
            if (t_num3==1)
            model_mesh_nei.nei_list_copoint.all{ii,1}=jj;
            else
                model_mesh_nei.nei_list_copoint.all{ii,1}=[model_mesh_nei.nei_list_copoint.all{ii,1} jj];
            end
            nei_pts_copoint_t{ii,t_num3}=temp1;
        end
%         if (t_num==4)
%             break; %no need search, save time;
%         end
    end
    for t1=1:t_num
        model_mesh_nei_tmp.nei_pts_coplane{ii,t1}=nei_pts_coplane{ii,t1};
        model_mesh_nei_tmp.nei_pts_coline{ii,t1}=nei_pts_coplane{ii,t1};
        model_mesh_nei_tmp.nei_pts_copoint{ii,t1}=nei_pts_coplane{ii,t1};
    end
    
    for t2=1:t_num2
        model_mesh_nei_tmp.nei_pts_coline{ii,t2+t_num}=nei_pts_coline_t{ii,t2};
        model_mesh_nei_tmp.nei_pts_copoint{ii,t2+t_num}=nei_pts_coline_t{ii,t2};
    end
    
    for t3=1:t_num3
        model_mesh_nei_tmp.nei_pts_copoint{ii,t3+t_num+t_num2}=nei_pts_copoint_t{ii,t3};
    end
    
    model_mesh_nei.nei_list_copoint.all{ii,1}=[model_mesh_nei.nei_list_coplane.all{ii,1} model_mesh_nei.nei_list_coline.all{ii,1} model_mesh_nei.nei_list_copoint.all{ii,1}];
    model_mesh_nei.nei_list_coline.all{ii,1}=[model_mesh_nei.nei_list_coplane.all{ii,1} model_mesh_nei.nei_list_coline.all{ii,1}];
    model_mesh_nei.tot_neighbour(ii,1)=t_num; %co-plane neighbours, 4 inner, 3 surface,  2 intercept, 1 corner
    model_mesh_nei.tot_neighbour(ii,2)=t_num+t_num2;%co-line neighbours
    model_mesh_nei.tot_neighbour(ii,3)=t_num+t_num2+t_num3;%co-point neighbours
    waitbar(ii/l, h);
end

for ii=1:l
    %get faces share with other tetrahedron
    pa=model_mesh_nei.t_sort(1,ii);
    pb=model_mesh_nei.t_sort(2,ii);
    pc=model_mesh_nei.t_sort(3,ii);
    pd=model_mesh_nei.t_sort(4,ii);    
    face{1,1}=[pa pb pc];face{1,2}=1;
    face{2,1}=[pa pb pd];face{2,2}=1;
    face{3,1}=[pa pc pd];face{3,2}=1;
    face{4,1}=[pb pc pd];face{4,2}=1;
    model_mesh_nei.nei_list_coplane.four_planes{ii,1}=face{1,1};
    model_mesh_nei.nei_list_coplane.four_planes{ii,2}=face{2,1};
    model_mesh_nei.nei_list_coplane.four_planes{ii,3}=face{3,1};
    model_mesh_nei.nei_list_coplane.four_planes{ii,4}=face{4,1};
    for jj=1:model_mesh_nei.tot_neighbour(ii,1)
        tmpa=[];tmpb=[];tmpc=[];tmpd=[];
        temp1=nei_pts_coplane{ii,jj};
        if temp1(1)==1 tmpa=pa; end
        if temp1(2)==1 tmpb=pb; end
        if temp1(3)==1 tmpc=pc; end
        if temp1(4)==1 tmpd=pd; end
        nei_list_coplane.co_face{ii,jj}=[tmpa tmpb tmpc tmpd];
        if (model_mesh_nei.tot_neighbour(ii,1)<=3)
            for chk_face=1:4
                if (nei_list_coplane.co_face{ii,jj}==face{chk_face,1})
                    face{chk_face,2}=0;
                end
            end
        end
    end
    
    if (model_mesh_nei.tot_neighbour(ii,1)<=3)
        t_num=0;
        for chk_face=1:4
            if (face{chk_face,2}==1)
                t_num=t_num+1;
                nei_list_coplane.surface{ii,t_num}=face{chk_face,1};
            end
        end
    end
    
end

for ii=1:l
    plane_nei_num=model_mesh_nei.nei_list_coplane.all{ii,1};
    for jj=1:4
        plane_neilist=model_mesh_nei.nei_list_coplane.four_planes{ii,jj};
        t_surface=1;
        for num=1:model_mesh_nei.tot_neighbour(ii,1);
            plane_comp=nei_list_coplane.co_face{ii,num};
            if isequal(plane_neilist,plane_comp)
                model_mesh_nei.nei_list_coplane.four_neighbour(ii,jj)=plane_nei_num(num);
                t_surface=0;
                %break;
            end
        end
        if (t_surface==1)
            model_mesh_nei.nei_list_coplane.four_neighbour(ii,jj)=nan; % indicate it is a surface, no neighbour
        end
    end
end

%find co-point neighbour for each point
for ii=1:l
    %get tetrahdron number share with this tetrahedron corner
    pa=model_mesh_nei.t_sort(1,ii);
    pb=model_mesh_nei.t_sort(2,ii);
    pc=model_mesh_nei.t_sort(3,ii);
    pd=model_mesh_nei.t_sort(4,ii);
    model_mesh_nei.nei_list_copoint.point{ii,1}=[pa pb pc pd];
    Ma=[];Mb=[];Mc=[];Md=[];
    for jj=1:model_mesh_nei.tot_neighbour(ii,3) %(ii,3) total neighbour points
        temp1=model_mesh_nei_tmp.nei_pts_copoint{ii,jj};
        if temp1(1)==1 Ma=[Ma model_mesh_nei.nei_list_copoint.all{ii,1}(1,jj)]; end
        if temp1(2)==1 Mb=[Mb model_mesh_nei.nei_list_copoint.all{ii,1}(1,jj)]; end
        if temp1(3)==1 Mc=[Mc model_mesh_nei.nei_list_copoint.all{ii,1}(1,jj)]; end
        if temp1(4)==1 Md=[Md model_mesh_nei.nei_list_copoint.all{ii,1}(1,jj)]; end
    end
    model_mesh_nei.nei_list_copoint.seperate{ii,1}=Ma;
    model_mesh_nei.nei_list_copoint.seperate{ii,2}=Mb;
    model_mesh_nei.nei_list_copoint.seperate{ii,3}=Mc;
    model_mesh_nei.nei_list_copoint.seperate{ii,4}=Md;  
end


%find co-point neighbour for each line
for ii=1:l
    %get tetrahdron number share with this tetrahedron line
    pa=model_mesh_nei.t_sort(1,ii);
    pb=model_mesh_nei.t_sort(2,ii);
    pc=model_mesh_nei.t_sort(3,ii);
    pd=model_mesh_nei.t_sort(4,ii);
    model_mesh_nei.nei_list_coline.line{ii,1}=[pa pb];
    model_mesh_nei.nei_list_coline.line{ii,2}=[pa pc];
    model_mesh_nei.nei_list_coline.line{ii,3}=[pa pd];
    model_mesh_nei.nei_list_coline.line{ii,4}=[pb pc];
    model_mesh_nei.nei_list_coline.line{ii,5}=[pb pd];
    model_mesh_nei.nei_list_coline.line{ii,6}=[pc pd];
    Mab=[];Mac=[];Mad=[];Mbc=[];Mbd=[];Mcd=[];
    for jj=1:model_mesh_nei.tot_neighbour(ii,2) %(ii,2) total neighbour lines
        temp1=model_mesh_nei_tmp.nei_pts_coline{ii,jj};
        if (temp1(1)==1 && temp1(2)==1) Mab=[Mab model_mesh_nei.nei_list_coline.all{ii,1}(1,jj)]; end
        if (temp1(1)==1 && temp1(3)==1) Mac=[Mac model_mesh_nei.nei_list_coline.all{ii,1}(1,jj)]; end
        if (temp1(1)==1 && temp1(4)==1) Mad=[Mad model_mesh_nei.nei_list_coline.all{ii,1}(1,jj)]; end
        if (temp1(2)==1 && temp1(3)==1) Mbc=[Mbc model_mesh_nei.nei_list_coline.all{ii,1}(1,jj)]; end
        if (temp1(2)==1 && temp1(4)==1) Mbd=[Mbd model_mesh_nei.nei_list_coline.all{ii,1}(1,jj)]; end
        if (temp1(3)==1 && temp1(4)==1) Mcd=[Mcd model_mesh_nei.nei_list_coline.all{ii,1}(1,jj)]; end
    end
    model_mesh_nei.nei_list_coline.seperate{ii,1}=Mab;
    model_mesh_nei.nei_list_coline.seperate{ii,2}=Mac;
    model_mesh_nei.nei_list_coline.seperate{ii,3}=Mad;
    model_mesh_nei.nei_list_coline.seperate{ii,4}=Mbc;
    model_mesh_nei.nei_list_coline.seperate{ii,5}=Mbd;
    model_mesh_nei.nei_list_coline.seperate{ii,6}=Mcd;
    %waitbar(ii/l, h);
end

for ii=1:size(model_mesh_nei.t,2)
    for jj=1:4
    	num=model_mesh_nei.nei_list_copoint.point{ii}(jj);
        model_mesh_nei.p_nei_tetra{num,1}=sort([ii model_mesh_nei.nei_list_copoint.seperate{ii,jj}]);
    end
end


if (chk==2) %save to file
    save(model_mesh.filename, 'model_mesh_nei')
end
close(h)

else
    load(model_mesh.filename, 'model_mesh_nei')
    t_load=model_mesh_nei.t;
    t_comp=model_mesh.t;
    if (isequal(t_load,t_comp)==1)
        disp('Neighbour list file is loaded.')
    else
        disp('Error! Load data is not consistent with current input! Nothing will be returned.')
        model_mesh_nei=[];
    end
end


toc
end