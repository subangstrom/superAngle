function [ dist, p_intercept, n_nei ] = get_intercept_model(sample_start_p, sample_start_num, line_direction, model_sample, check )
%Weizong Xu, August, 2017

line_point=sample_start_p;
mesh_num=sample_start_num;
if (check>0) %in forward direction
    while (isnan(mesh_num)==0)
        n_nei=mesh_num;
        [  p_forward, n_forward, ~, ~, ~ ] = intercept_line_tetra_fast( line_point, line_direction, mesh_num, model_sample.p,model_sample.t_sort, model_sample.nei_list_coplane);
        %[  p_forward, n_forward, ~, ~, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_sample );
        line_point=p_forward;
        mesh_num=n_forward;
        %sample_list_all=[sample_list_all n_forward n_backward];
    end
    p_intercept=line_point;
    d_intercept=line_point-sample_start_p;
    %dist=norm(line_point-sample_start_p);
    dist=sqrt(d_intercept(1)^2+d_intercept(2)^2+d_intercept(3)^2);
end

if (check<0) %in backward direction
    line_point=sample_start_p;
    mesh_num=sample_start_num;
    while (isnan(mesh_num)==0)
        n_nei=mesh_num;
        [  ~, ~, p_backward, n_backward, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_sample );
        line_point=p_backward;
        mesh_num=n_backward;
        %sample_list_all=[sample_list_all n_forward n_backward];
    end
    %line_start=line_point;
    p_intercept=line_point;
    dist=norm(line_point-sample_start_p);
end

if (check==0) %full range search
    sample_list_all=[sample_start_num];
    line_point=sample_start_p;
    mesh_num=sample_start_num;
    while (isnan(mesh_num)==0)
        [  p_forward, n_forward, ~, n_backward, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_sample );
        line_point=p_forward;
        mesh_num=n_forward;
        sample_list_all=[sample_list_all n_forward n_backward];
    end
    p_intercept.p_end=line_point;

    line_point=sample_start_p;
    mesh_num=sample_start_num;
    while (isnan(mesh_num)==0)
        [  ~, n_forward, p_backward, n_backward, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_sample );
        line_point=p_backward;
        mesh_num=n_backward;
        sample_list_all=[sample_list_all n_forward n_backward];
    end
    p_intercept.p_start=line_point;
    dist=norm(p_intercept.p_end-p_intercept.p_start);
    sample_list_all(isnan(sample_list_all))=[];
    n_nei=unique(sample_list_all); %get neighbour list of sample that line is possibily intercept
end

end

