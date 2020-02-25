%Membrane constants that display upper and lower bounds

%Publication #1: Chemical Reacting Flow Theory and Simulations 
%                Author Robert J Kee
%Publication #2: Compilation of mechanical properties for the structural 
%                analysis of solid oxide fuel cell stacks. Constitutive 
%                materials of anode-supported cells, Arata Nakajo
%Publication #3: A modeling study on concentration overpotentials of a
%                reversible solid oxide fuel cell, Meng Ni


%Anode (LSM + Argon + Oxygen)
LSM_Porosity = mean([.04 .09 .12 .29 .30 .35]);
LSM_Torosity = mean([4 6]);
LSM_radius = .625; %um

%Cathode (NIYSZ + Hydrogen + Water)
NIYSZ_Porosity  = mean([.26 .35 .40]);
NIYSZ_Torosity =  mean([4.5 6]);
NIYSZ_Spec_TPB = 2.40*10^9; % 1/cm^2
NIYSZ_radius = .5; %um
