function [delta_O2_eff] = O_in_Ar_delta(To, Po, L_P, L_T, L_R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
LSM_radius = L_R;
LSM_Torosity = L_T;
LSM_Porosity = L_P;

%Publication #1: A modeling study on concentration overpotentials of a
%reversible solid oxide fuel cell Meng Ni, 

R = 8.314459;  %J?mol?1?K?1
kb = 1.38066*10^-23; %J/K
M_O2 = 0.016; %kg/mol

%-----------------------------------------------------------------
Dp = 9.77*10^(-6)*To^1.736;
%constants use pressure in atm 
%hydrogen diffusion is supercritical water requires a mixed density while
%oxygen in argon requires pressure of the system

%1 MPa = 9.869 atm
Patm = Po*9.869;
O2_D  = Dp/Patm;
%cm^2/s diffusion unit
%-----------------------------------------------------------------
D_o2_k  = (10^(-2))*(4/3)*LSM_radius*sqrt(8*R*To/(3.14159*M_O2));

delta_O2_eff = D_o2_k/(O2_D + D_o2_k);
end