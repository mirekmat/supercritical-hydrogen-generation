function [n_potential] = nconH2_function(J, To, Po, Rxn, thickness_NIYSZ, ...
                                   NIYSZ_Porosity, NIYSZ_Torosity, NIYSZ_radius)
%Concentration Potential Anode 
%   Detailed explanation goes here

t = thickness_NIYSZ;
D_HW = W_in_H_diffusion(To, Po, NIYSZ_Porosity, NIYSZ_Torosity, NIYSZ_radius);

F = 96485.33289; %C mol?1
R = 8.314459;  %J?mol?1?K?1
%Pressure changes 1 MPA to J/cm^3
% To  Units [K]
% J   Units [A/cm^2];
% Rxn Units [0-1]
var = (R*To*J*t)/(2*F*D_HW);

x1 = 1 + var/(Po*10^6*Rxn);
x2 = 1 - var/(Po*10^6*(1-Rxn));

x = x1/x2;

   if (x < 0)
      n_potential = NaN;
   else 
      n_potential = (R*To/(2*F))*reallog(x);
   end
end