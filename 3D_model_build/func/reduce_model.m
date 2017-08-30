function [ model_out ] = reduce_model( model_in )
%reduce model size for storage, only save data that necessary for
%post calculation
%Weizong, August, 2017

model_out.p=model_in.p;
model_out.t=model_in.t;
model_out.t_sort=model_in.t_sort;
model_out.nei_list_coplane.four_neighbour=model_in.nei_list_coplane.four_neighbour;
model_out.tetra_volume=model_in.tetra_volume;
model_out.start_point=model_in.start_point;
model_out.start_tetra_num=model_in.start_tetra_num;
model_out.filename=model_in.filename;

end

