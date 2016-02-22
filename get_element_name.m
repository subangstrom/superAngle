function [ sym_out, sym, sym_shell] = get_element_name( element_number, element_shell )
%Get element name
%Weizong Xu, wxu4@ncsu.edu, April, 2015
switch element_number
    case 1
        sym='H';
    case 2
        sym='He';
    case 3
        sym='Li';
    case 4
        sym='Be';
    case 5
        sym='B';
    case 6
        sym='C';
    case 7
        sym='N';
    case 8
        sym='O';
    case 9
        sym='F';
    case 10
        sym='Ne';
    case 11
        sym='Na';
    case 12
        sym='Mg';
    case 13
        sym='Al';
    case 14
        sym='Si';
    case 15
        sym='P';
    case 16
        sym='S';
    case 17
        sym='Cl';
    case 18
        sym='Ar';
    case 19
        sym='K';
    case 20
        sym='Ca';
    case 21
        sym='Sc';
    case 22
        sym='Ti';
    case 23
        sym='V';
    case 24
        sym='Cr';
    case 25
        sym='Mn';
    case 26
        sym='Fe';
    case 27
        sym='Co';
    case 28
        sym='Ni';
    case 29
        sym='Cu';
    case 30
        sym='Zn';
    case 31
        sym='Ga';
    case 32
        sym='Ge';
    case 33
        sym='As';
    case 34
        sym='Se';
    case 35
        sym='Br';
    case 36
        sym='Kr';
    case 37
        sym='Rb';
    case 38
        sym='Sr';
    case 39
        sym='Y';
    case 40
        sym='Zr';
    case 41
        sym='Nb';
    case 42
        sym='Mo';
    case 43
        sym='Tc';
    case 44
        sym='Ru';
    case 45
        sym='Rh';
    case 46
        sym='Pd';
    case 47
        sym='Ag';
    case 48
        sym='Cd';
    case 49
        sym='In';
    case 50
        sym='Sn';
    case 51
        sym='Sb';
    case 52
        sym='Te';
    case 53
        sym='I';
    case 54
        sym='Xe';
    case 55
        sym='Cs';
    case 56
        sym='Ba';
    case 57
        sym='La';
    case 58
        sym='Ce';
    case 59
        sym='Pr';
    case 60
        sym='Nd';
    case 61
        sym='Pm';
    case 62
        sym='Sm';
    case 63
        sym='Eu';
    case 64
        sym='Gd';
    case 65
        sym='Tb';
    case 66
        sym='Dy';
    case 67
        sym='Ho';
    case 68
        sym='Er';
    case 69
        sym='Tm';
    case 70
        sym='Yb';
    case 71
        sym='Lu';
    case 72
        sym='Hf';
    case 73
        sym='Ta';
    case 74
        sym='W';
    case 75
        sym='Re';
    case 76
        sym='Os';
    case 77
        sym='Ir';
    case 78
        sym='Pt';
    case 79
        sym='Au';
    case 80
        sym='Hg';
    case 81
        sym='Tl';
    case 82
        sym='Pb';
    case 83
        sym='Bi';
    case 84
        sym='Po';
    case 85
        sym='At';
    case 86
        sym='Rn';
    case 87
        sym='Fr';
    case 88
        sym='Ra';
    case 89
        sym='Ac';
    case 90
        sym='Th';
    case 91
        sym='Pa';
    case 92
        sym='U';
    case 93
        sym='Np';
    case 94
        sym='Pu';
    otherwise
        sym='?';
end

switch element_shell
    case 1
        sym_shell='K';
    case 2
        sym_shell='L';
    case 3
        sym_shell='M';
    case 4
        sym_shell='N';
    case 5
        sym_shell='O';
    otherwise
        sym='?';
end

sym_out = [sym,'-',sym_shell];

end

