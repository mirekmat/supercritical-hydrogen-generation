%% Specific Pressure and Temperature

    temp_indiv = 273+700;  %.C
    press_indiv = 400; %MPa

    B_Dissociation_indiv = 0;
    B_Ionic_indiv        = 0;
    KINETIC_P_indiv      = 0;
    P_O2c_indiv          = 0;
    J_indiv              = 0;
    NONIDEAL_P_indiv     = 0;
    IDEAL_P_indiv        = 0;
    nohmic_a             = 0;
    nact_a               = 0;
    nconH2_a             = 0;
    nconO2_a             = 0;
    
    rxn = .5;
    %Setting Variables
    %===============================================================================
    row = round((temp_indiv - 20 - 273.15)/5);
    col = round(press_indiv - 22);
    P_indiv              = POTENTIAL(col, row);
    %Calculating K Factor
    %===============================================================================
    K_indiv = exp(GIBBS_FORMATION(col,row)*1000/(8.314*temp_indiv));
    
    %Electronic Region Check
    %===============================================================================
    Min_950 = 10^(-16);
    Max_950 = 10^7;
    Min_650 = 10^(-31);
    Max_650 = 10^12;
      
    Min_l = Min_650 + ((Min_950 - Min_650)/(950-650))*(temp_indiv - (650 + 273.15));
    Max_l = Max_650 + ((Max_950 - Max_650)/(950-650))*(temp_indiv - (650 + 273.15));
    
    E_Dissoc_indiv = -2.84215 + 5.02487*10^(-4)*temp_indiv -1.23698*10^(-8)*temp_indiv^2;
    
    for J_index = -5:.01:5
      
      %Mapping Thermodynamic Potential
      %===============================================================================
      [pot_ideal, non_ideal] = thermo__potentialb(rxn,OXYGEN_ANODE,temp_indiv,press_indiv,K_indiv);
      IDEAL_P_indiv   (J_index) = pot_ideal;
      NONIDEAL_P_indiv(J_index) = non_ideal;
      
      %Electronic Region Check
      %===============================================================================
      P_O2c_indiv(J_index) = sqrt(K_indiv)*(1 - rxn)/rxn;
      
      if(P_O2c_indiv(J_index) < Min_l || P_O2c_indiv(J_index) > Max_l)
          B_Ionic_indiv (J_index) = 0;
      else
          B_Ionic_indiv (J_index) = 1;
      end
      
      %Zirconia Dissociation Check
      %===============================================================================
      
      if(E_Dissoc_indiv < NONIDEAL_P_indiv(J_index))
          B_Dissociation_indiv (J_index) = 1;
      else
          B_Dissociation_indiv (J_index) = 0;
      end
      
      J_indiv(J_index) = 3.05; %A/cm2 current density intitial convergence value
      Pre_kinetic = NONIDEAL_P_indiv(J_index);
      KINETIC_P_indiv(J_index) = -100; %just putting a case that would initiate
      %                               the while loop. 
      %Mapping Kinetic Potential
      %===============================================================================
      [nohmic, nact, nconH2, nconO2, kinetic_pot] = kinetics(J_indiv(J_index), temp_indiv, press_indiv, rxn, ...
                               YSZ_T*10000, NIYSZ_T, LSM_T, ...
                               NIYSZ_Porosity, NIYSZ_Torosity, NIYSZ_radius, ...
                               LSM_Porosity, LSM_Torosity, LSM_radius, OXYGEN_ANODE);
      KINETIC_P_indiv(J_index) = kinetic_pot;
      nohmic_a(J_index)  = nohmic;
      nact_a(J_index)    = nact;
      nconH2_a(J_index)  = nconH2;
      nconO2_a(J_index)  = nconO2;
                           
      end   
    end
    
%% Illustration of current + voltage subcomponents 
figure(20)
X_rxn = 1/10:1/10:999/10;

Pot       = zeros(999,1);
Zir       = zeros(999,1);
for i = 1:length(Pot) 
    % set each element to 0
    Pot(i) = P_indiv;
    Zir(i) = E_Dissoc_indiv;
end

yyaxis left
title(['Voltage + Current at ' , num2str(temp_indiv), ' K & ' , num2str(press_indiv) , ' MPa'])
plot(X_rxn, Pot,'-g');
hold on
plot(X_rxn, NONIDEAL_P_indiv, '-r');
hold on
plot(X_rxn, NONIDEAL_P_indiv + KINETIC_P_indiv, '-b');
hold on
plot(X_rxn, NONIDEAL_P_indiv - nohmic_a, '--r');
hold on
plot(X_rxn, NONIDEAL_P_indiv - nohmic_a - nact_a, '--b');
hold on
plot(X_rxn, NONIDEAL_P_indiv - nohmic_a - nact_a - nconH2_a, '--m');
hold on
plot(X_rxn, NONIDEAL_P_indiv - nohmic_a - nact_a - nconH2_a - nconO2_a, '--c');
hold on
l2 = plot(X_rxn,Zir, ':r');
set(l2,'linewidth',2);
hold on

xlabel('Reaction Completion [%]') % x-axis label
ylabel('Potential [V]') % y-axis label
set(gca,'XLim',[0 100])
set(gca,'XTick',(0:10:100))
set(gca,'YLim',[-2.5 0])
set(gca,'YTick',((-2.5):.5:0))

yyaxis right
l1 = plot(X_rxn,J_indiv,'-b');
set(l1,'linewidth',2);
ylabel('Current [A/cm^2]') % y-axis label
set(gca,'YLim',[0 4])
set(gca,'YTick',(0:.5:4))
legend('Thermodynamic','Non-Ideal Thermo Addition', 'Kinetic Addition',...
'Ohmic Resistance', 'Activation Polarization', 'Hydrogen Concentration Polar.',...
'Oxygen Concentration Polar.','Zirconia Dissociation','Current Density')
hold off

%% I-V CURVES

    KINETIC_P_indiv2      = zeros(100,3);
    P_O2c_indiv2          = zeros(100,3);
    J_indiv2              = zeros(100,3);
    NONIDEAL_P_indiv2     = zeros(100,3);
    IDEAL_P_indiv2        = zeros(100,3);
    nohmic_a2             = zeros(100,3);
    nact_a2               = zeros(100,3);
    nconH2_a2             = zeros(100,3);
    nconO2_a2             = zeros(100,3);
    
    
for J_index = 1:1:100
    
    for rxn_index2 = 1:1:3
           
           J2   = (J_index-1)/10;
           rxn2 = ((rxn_index2-1)*49 + 1)/100;
           
           [pot_ideal, non_ideal] = thermo__potentialb(rxn2,OXYGEN_ANODE,temp_indiv,press_indiv,K_indiv);
           IDEAL_P_indiv2   (J_index, rxn_index2) = pot_ideal;
           NONIDEAL_P_indiv2(J_index, rxn_index2) = non_ideal;
      
     
           %Mapping Kinetic Potential
           %===============================================================================
           [nohmic, nact, nconH2, nconO2, kinetic_pot] = kinetics(J2, temp_indiv, press_indiv, rxn2, ...
                               YSZ_T*10000, NIYSZ_T, LSM_T, ...
                               NIYSZ_Porosity, NIYSZ_Torosity, NIYSZ_radius, ...
                               LSM_Porosity, LSM_Torosity, LSM_radius, OXYGEN_ANODE);
                           
            KINETIC_P_indiv2(J_index, rxn_index2) = kinetic_pot;
            nohmic_a2(J_index, rxn_index2)  = nohmic;
            nact_a2(J_index, rxn_index2)  = nact;
            nconH2_a2(J_index, rxn_index2)  = nconH2;
            nconO2_a2(J_index, rxn_index2)  = nconO2; 
            
   
end