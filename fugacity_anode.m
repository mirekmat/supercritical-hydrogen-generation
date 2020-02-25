function [Ar_theta,O2_theta, Z] = fugacity_anode(Ar_y,O2_y,To, Po)
%Fugacity_anode.m 
% Thermodynamic Modeling of the Water-Gas Shift Reaction in 
% Supercritical Water for Hydrogen Production


%Ar Fugacity Data
Ar_Tc = 151;                          %K
Ar_Pc = 4.86;                         %MPa
Ar_Vc = .291*151*8.314/(4.86*10^6);   %m3/mol
Ar_w  = 0;                            %w

%O2 Fugacity Data
O2_Tc = 155;
O2_Pc = 5.08;
O2_Vc = 0.288*155*8.314/(5.08*10^6);
O2_w  = 0.021;

alpha_O2  = (1+(0.37464 + 1.5422*O2_w - 0.26992*O2_w^2)*(1-sqrt(To/O2_Tc)))^2;
alpha_Ar  = (1+(0.37464 + 1.5422*Ar_w - 0.26992*Ar_w^2)*(1-sqrt(To/Ar_Tc)))^2;

A_O2  = 0.45726*(8.314^2)*(O2_Tc)^2*alpha_O2/(O2_Pc*10^6);
A_Ar  = 0.45726*(8.314^2)*(Ar_Tc)^2*alpha_Ar/(Ar_Pc*10^6);

B_O2  = 0.0778*8.3154*O2_Tc/(O2_Pc*10^6);
B_Ar  = 0.0778*8.3154*Ar_Tc/(Ar_Pc*10^6);

B_m =  B_O2*O2_y + B_Ar*Ar_y;

k_oa = 1 - 8*(Ar_Vc*O2_Vc)^2/(((Ar_Vc)^(1/3))*((O2_Vc)^(1/3)))^3; 
k_oo = 1 - 8*(O2_Vc*O2_Vc)^2/(((O2_Vc)^(1/3))*((O2_Vc)^(1/3)))^3; 
k_aa = 1 - 8*(Ar_Vc*Ar_Vc)^2/(((Ar_Vc)^(1/3))*((Ar_Vc)^(1/3)))^3; 

A_oa = (1-k_oa)*sqrt(A_Ar*A_O2);
A_oo = (1-k_oo)*sqrt(A_O2*A_O2);
A_aa = (1-k_aa)*sqrt(A_Ar*A_Ar);

A_ao = A_oa;

A_m = A_aa*Ar_y^2  + A_ao*O2_y*Ar_y ...
    + A_oo*O2_y^2  + A_ao*O2_y*Ar_y;

A =  A_m*Po*10^6/(8.314^2*To^2);
B =  B_m*Po*10^6/(8.314*To);

p = [1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)];

r = roots(p);
Z = max(r);

Ar_q  = 2*(O2_y*A_ao + Ar_y*A_aa)/A_m - B_Ar/B;
O2_q  = 2*(O2_y*A_oo + Ar_y*A_oa)/A_m - B_O2/B;

Term_z = reallog((Z + (1+sqrt(2)*B))/(Z - (1-sqrt(2)*B)));

Ln_Ar_theta  = (B_Ar/B_m)*(Z-1)  - reallog(Z-B) - (A/(2*B*sqrt(2)))*Ar_q*Term_z;
Ln_O2_theta  = (B_O2/B_m)*(Z-1)  - reallog(Z-B) - (A/(2*B*sqrt(2)))*O2_q*Term_z;

Ar_theta  = exp(Ln_Ar_theta);
O2_theta  = exp(Ln_O2_theta);

end

