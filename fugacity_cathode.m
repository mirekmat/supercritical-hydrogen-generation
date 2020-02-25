function [H2O_theta,H2_theta,O2_theta, Z] = fugacity_cathode(H2O_y,H2_y,O2_y,To, Po)
%Fugacity_cathode.m 
% Thermodynamic Modeling of the Water-Gas Shift Reaction in 
% Supercritical Water for Hydrogen Production


%H2 Fugacity Data
H2_Tc = 33.19;        %K
H2_Pc = 1.31;         %MPa
H2_Vc = 6.42*10^-5;   %m3/mol
H2_w  = -.2320;       %w

%H2O Fugacity Data
H2O_Tc = 647.13;
H2O_Pc = 22.06;
H2O_Vc = 5.6*10^-5;
H2O_w  = .3449;

%O2 Fugacity Data
O2_Tc = 155;
O2_Pc = 5.08;
O2_Vc = 0.288*155*8.314/(5.08*10^6);
O2_w  = 0.021;

alpha_O2  = (1+(0.37464 + 1.5422*O2_w - 0.26992*O2_w^2)*(1-sqrt(To/O2_Tc)))^2;
alpha_H2O = (1+(0.37464 + 1.5422*H2O_w - 0.26992*H2O_w^2)*(1-sqrt(To/H2O_Tc)))^2;
alpha_H2  = (1+(0.37464 + 1.5422*H2_w - 0.26992*H2_w^2)*(1-sqrt(To/H2_Tc)))^2;

A_H2O = 0.45726*(8.314^2)*(H2O_Tc)^2*alpha_H2O/(H2O_Pc*10^6);
A_O2  = 0.45726*(8.314^2)*(O2_Tc)^2*alpha_O2/(O2_Pc*10^6);
A_H2  = 0.45726*(8.314^2)*(H2_Tc)^2*alpha_H2/(H2_Pc*10^6);

B_H2O = 0.0778*8.3154*H2O_Tc/(H2O_Pc*10^6);
B_O2  = 0.0778*8.3154*O2_Tc/(O2_Pc*10^6);
B_H2  = 0.0778*8.3154*H2_Tc/(H2_Pc*10^6);

B_m = B_H2O*H2O_y + B_O2*O2_y + B_H2*H2_y;

k_wh = 1 - 8*(H2O_Vc*H2_Vc)^2/(((H2O_Vc)^(1/3))*((H2_Vc)^(1/3)))^3;
k_wo = 1 - 8*(H2O_Vc*O2_Vc)^2/(((H2O_Vc)^(1/3))*((O2_Vc)^(1/3)))^3; 
k_ho = 1 - 8*(H2_Vc*O2_Vc)^2/(((H2_Vc)^(1/3))*((O2_Vc)^(1/3)))^3; 
k_ww = 1 - 8*(H2O_Vc*H2O_Vc)^2/(((H2O_Vc)^(1/3))*((H2O_Vc)^(1/3)))^3;
k_oo = 1 - 8*(O2_Vc*O2_Vc)^2/(((O2_Vc)^(1/3))*((O2_Vc)^(1/3)))^3; 
k_hh = 1 - 8*(H2_Vc*H2_Vc)^2/(((H2_Vc)^(1/3))*((H2_Vc)^(1/3)))^3; 

A_wh = (1-k_wh)*sqrt(A_H2O*A_H2);
A_wo = (1-k_wo)*sqrt(A_H2O*A_O2);
A_ho = (1-k_ho)*sqrt(A_H2*A_O2);
A_ww = (1-k_ww)*sqrt(A_H2O*A_H2O);
A_oo = (1-k_oo)*sqrt(A_O2*A_O2);
A_hh = (1-k_hh)*sqrt(A_H2*A_H2);

A_hw = A_wh;
A_ow = A_wo;
A_oh = A_ho;

A_m = A_ww*H2O_y^2 + A_wh*H2O_y*H2_y + A_wo*O2_y*H2O_y ...
    + A_hh*H2_y^2  + A_ho*O2_y*H2_y  + A_wh*H2O_y*H2_y ...
    + A_oo*O2_y^2  + A_ho*O2_y*H2_y  + A_wo*O2_y*H2O_y;


A =  A_m*Po*10^6/(8.314^2*To^2);
B =  B_m*Po*10^6/(8.314*To);

p = [1 -(1-B) (A-2*B-3*B^2) -(A*B-B^2-B^3)];

r = roots(p);
Z = max(r);

H2O_q = 2*(H2O_y*A_ww + O2_y*A_wo + H2_y*A_wh)/A_m - B_H2O/B;
H2_q  = 2*(H2O_y*A_hw + O2_y*A_ho + H2_y*A_hh)/A_m - B_H2/B;
O2_q  = 2*(H2O_y*A_ow + O2_y*A_oo + H2_y*A_oh)/A_m - B_O2/B;

Term_z = reallog((Z + (1+sqrt(2)*B))/(Z - (1-sqrt(2)*B)));

Ln_H2O_theta = (B_H2O/B_m)*(Z-1) - reallog(Z-B) - (A/(2*B*sqrt(2)))*H2O_q*Term_z;
Ln_H2_theta  = (B_H2/B_m)*(Z-1)  - reallog(Z-B) - (A/(2*B*sqrt(2)))*H2_q*Term_z;
Ln_O2_theta  = (B_O2/B_m)*(Z-1)  - reallog(Z-B) - (A/(2*B*sqrt(2)))*O2_q*Term_z;

H2O_theta = exp(Ln_H2O_theta);
H2_theta  = exp(Ln_H2_theta);
O2_theta  = exp(Ln_O2_theta);

end

