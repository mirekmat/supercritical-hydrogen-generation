function [nohmic, nact, nconH2, nconO2, kinetic_pot] = kinetics(J, To, Po, Rxn, ...
                         YSZ_t, NIYSZ_t, LSM_t,   ...
                         NIYSZ_Porosity, NIYSZ_Torosity, NIYSZ_radius, ...
                         LSM_Porosity, LSM_Torosity, LSM_radius, o_ratio)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

nohmic = nohmic_function(J, YSZ_t, To); %thickness (um)

nact   = nact_function(To, J, Rxn); 

nconH2 = nconH2_function(J, To, Po, Rxn, NIYSZ_t, ...
                            NIYSZ_Porosity, NIYSZ_Torosity, NIYSZ_radius);
                        
nconO2 = nconO2_function(J, To, Po, o_ratio, LSM_t, ...
                            LSM_Porosity, LSM_Torosity, LSM_radius); 

kinetic_pot = (nohmic + nact + nconH2 + nconO2); 

end

