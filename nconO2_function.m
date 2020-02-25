function [n_potential] = nconO2_function(J, To, Po, o_ratio, thickness_LSM, ...
                                 LSM_Porosity, LSM_Torosity, LSM_radius)
                             
F = 96485.33289; %C mol?1
R = 8.314459;  %J?mol?1?K?1
Pc = Po; % pressure cathode MPa
Pobulk = Po*o_ratio;


L_P = LSM_Porosity;
L_T = LSM_Torosity;
L_R = LSM_radius;

delta_oxygen = O_in_Ar_delta(To, Po, L_P, L_T, L_R);
diffusion_oxygen = O_in_Ar_diffusion(To, Po, L_P, L_T, L_R);

exponential_term = exp((R*To/(4*F))*(-delta_oxygen*thickness_LSM/(diffusion_oxygen*Pc))*J);
top = Pc/delta_oxygen - ((Pc/delta_oxygen)-Pobulk)*exponential_term;
              
x = top/Pobulk;

   if (x < 0)
      n_potential = 100;
   else 
      n_potential = (R*To/(4*F))*reallog(x);
   end                                                     
end