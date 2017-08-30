function [ shift_vector] = find_center( atomsP )
%find center of the lattice
%@@@@record center point data for unit cell cut@@@@
%Design by Weizong 2011 in FORTRAN, rewritten to Matlab m format.

tout1=size(atomsP);

unit1=[1e20,1e20,1e20];
unit2=[0,0,0];


for i=1:tout1
    if (atomsP(i,1)>unit2(1)) 
        unit2(1)=atomsP(i,1);
    end
    
    if (atomsP(i,2)>unit2(2)) 
        unit2(2)=atomsP(i,2);
    end
    
    if (atomsP(i,3)>unit2(3))
        unit2(3)=atomsP(i,3);
    end

    if (atomsP(i,1)<unit1(1)) 
        unit1(1)=atomsP(i,1);
    end
    
    if (atomsP(i,2)<unit1(2))
        unit1(2)=atomsP(i,2);
    end
    
    if (atomsP(i,3)<unit1(3))
        unit1(3)=atomsP(i,3);
    end
end

disp('Running... Finding the center of atom model')
disp(['min point ', num2str(unit1)]);
disp(['max point ', num2str(unit2)]);

unit2=0.5*(unit1+unit2);   % find the middle

%initiate value
i_min=1;
distance_min=(atomsP(1,1)-unit2(1))^2+(atomsP(1,2)-unit2(2))^2+(atomsP(1,3)-unit2(3))^2;

for i=1:tout1
    
    distance=(atomsP(i,1)-unit2(1))^2+(atomsP(i,2)-unit2(2))^2+(atomsP(i,3)-unit2(3))^2;
    if (distance<distance_min)
        distance_min=distance;
        i_min=i;
    end
end

disp(['Find a point! Distance to center:',num2str(sqrt(distance_min))]); 
disp([' @ ', ' n= ',num2str(i_min), ' x= ', num2str(atomsP(i_min,1)),' y= ',num2str(atomsP(i_min,2)),' z= ',num2str(atomsP(i_min,3))]);


shift_vector(1)=atomsP(i_min,1);
shift_vector(2)=atomsP(i_min,2);
shift_vector(3)=atomsP(i_min,3);
disp ('Done. Return center coordinate')
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


end

