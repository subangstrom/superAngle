function [ atomsP_Sr_out, atomsP_Ti_out, atomsP_O_out ] = lattice_translation( atomsP_Sr, atomsP_Ti, atomsP_O, shift_vector )
%Weizong Xu, August, 2017

atomsP_Sr_out(:,1)=atomsP_Sr(:,1)+shift_vector(1);
atomsP_Sr_out(:,2)=atomsP_Sr(:,2)+shift_vector(2);
atomsP_Sr_out(:,3)=atomsP_Sr(:,3)+shift_vector(3);

atomsP_Ti_out(:,1)=atomsP_Ti(:,1)+shift_vector(1);
atomsP_Ti_out(:,2)=atomsP_Ti(:,2)+shift_vector(2);
atomsP_Ti_out(:,3)=atomsP_Ti(:,3)+shift_vector(3);

atomsP_O_out(:,1)=atomsP_O(:,1)+shift_vector(1);
atomsP_O_out(:,2)=atomsP_O(:,2)+shift_vector(2);
atomsP_O_out(:,3)=atomsP_O(:,3)+shift_vector(3);

end

