function [D_H2O_eff] = W_in_H_diffusion(To, Po, N_P, N_T, N_R)

NIYSZ_radius = N_R;
NIYSZ_Torosity = N_T;
NIYSZ_Porosity = N_P;

%Publication #1: A modeling study on concentration overpotentials of a
%reversible solid oxide fuel cell Meng Ni, 

R = 8.314459;  %J?mol?1?K?1
kb = 1.38066*10^-23; %J/K
M_H2O = 0.018; %kg/mol
M_H2  = 0.002; %kg/mol

Char_L_H2  = 2.827; %A
Char_L_H2O = 2.641; %A
Char_L = .5*(Char_L_H2O + Char_L_H2);

Lenard_Jones_H2  = 5.144; %meV
Lenard_Jones_H2O = 69.72; %meV
LJE = (Lenard_Jones_H2O*Lenard_Jones_H2)^.5/(6.242*10^12); %J
tho = kb*To/LJE;

CI_A0 = 1.06036;
CI_A1 = 0.15610;
CI_A2 = 0.19300;
CI_A3 = 0.47635;
CI_A4 = 1.03587;
CI_A5 = 1.52996;
CI_A6 = 1.76474;
CI_A7 = 3.89411;

omega = CI_A0/tho^CI_A1 + CI_A2/exp(CI_A3*tho) + CI_A4/exp(CI_A5*tho) +...
        CI_A6/(CI_A7*tho);
    
M_red    = sqrt((1000/M_H2O + 1000/M_H2));

D_h2o_h2 = 0.001858*M_red*To^1.5/(omega*Char_L^2*Po*9.86923);

D_h2o_k  = (10^(-4))*(4/3)*NIYSZ_radius*sqrt(8*R*To/(3.14159*M_H2O));

iD_H2O_eff = (NIYSZ_Torosity/NIYSZ_Porosity)*(1/D_h2o_h2 + 1/D_h2o_k);

D_H2O_eff = 1/iD_H2O_eff;

end

