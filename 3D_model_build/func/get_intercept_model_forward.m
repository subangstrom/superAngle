function [ dist ] = get_intercept_model_forward(sample_start_p, sample_start_num, line_direction, p, t_sort, four_neighbour)
%Weizong Xu, August, 2017

line_point=sample_start_p;
mesh_num=sample_start_num;
%if (check>0) %in forward direction
    while (isnan(mesh_num)==0)
        [  p_forward, n_forward ] = intercept_line_tetra_fast_forward( line_point, line_direction, mesh_num, p, t_sort, four_neighbour);
        %[  p_forward1, n_forward1, ~, ~, ~] = intercept_line_tetra_fast( line_point, line_direction, mesh_num, model_sample.p,model_sample.t_sort, model_sample.nei_list_coplane);
        %[  p_forward, n_forward, ~, ~, ~ ] = intercept_line_tetra( line_point, line_direction, mesh_num, model_sample );
        %if (isequal(p_forward,p_forward1)==0)
        %    disp('Something is wrong!');
        %    [  p_forward, n_forward, ~] = intercept_line_tetra_fast_forward( line_point, line_direction, mesh_num, model_sample.p,model_sample.t_sort, model_sample.nei_list_coplane);
        %end
        line_point=p_forward;
        mesh_num=n_forward;
    end
    d_intercept=line_point-sample_start_p;
    %dist=norm(line_point-sample_start_p);
    dist=sqrt(d_intercept(1)^2+d_intercept(2)^2+d_intercept(3)^2);
%end



end

