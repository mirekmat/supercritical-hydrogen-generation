function [H2_D] = diffusion_hydrogen(T,p)

%water density should be in kg/m^3
%diffusion constant m^2/s

%In this study, the SPC/E model was used to perform classical
%molecular dynamics simulations of pure water as well as
%infinitely dilute aqueous solutions of O2, H2, and an OH radical
%species, in the supercritical region, as well as along the liquid
%branch of the liquid?vapor coexistence line. The transport
%properties of pure water and these infinitely dilute species were
%calculated using the density expansion of the correlation
%function expression for the diffusion coefficients of threedimensional
%gases, based on a hard sphere collision model, as
%formulated by Kawaski and Oppeinheim.

H2_a  = 116.211;
H2_alpha = 0.0538917;
H2_b1 = 1.0;
H2_b2 = 372312;
H2_b3 = -1794.21;
H2_b4 = 1.42193;
H2_c1 = 1.0;
H2_c2 = 476427;
H2_c3 = -1880.99;
H2_c4 = 1.62233;
H2_d1 = 1.0;
H2_d2 = -377905;
H2_d3 = 1654.15;
H2_d4 = -1.39764;
 
Dp = H2_a*T^H2_alpha + p(H2_b1*T^(-2) + H2_b2*T^(-1) + H2_b3 + H2_b4*T) + ...
     p^2*reallog(p)*(H2_c1*T^(-2) + H2_c2*T^(-1) + H2_c3 + H2_c4*T) + ...
     p^2*(H2_d1*T^(-2) + H2_d2*T^(-1) + H2_d3 + H2_d4*T);

H2_D  = Dp/p;

end



