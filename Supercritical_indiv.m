
function [] = Supercritical_indiv(temp_indiv, press_indiv, LSM_T, YSZ_T, NIYSZ_T, OXYGEN_ANODE)
%% Load Thermo Data + Membrane Data

load('C:\Users\Matthew\Documents\Supercritical Water\Thermodynamics\thermodynamic_data.mat')
GIBBS_FORMATION = -2*W_GIBBS + (O_GIBBS + 2*H_GIBBS)+50.30007462*2 +237.141*2;
POTENTIAL = -GIBBS_FORMATION/(4*96.485);
membrane_characteristics

filename = "T" + string(temp_indiv) + "P" + string(press_indiv) + ...
           "L" + string(LSM_T) + "Y" + string(YSZ_T) + "N" + string(NIYSZ_T) + ...
           "O" + string(OXYGEN_ANODE);
       
filename = strrep(filename,'.','_');
       
%% Specific Pressure and Temperature
    
    %%Global Variables

    %Structural/Geometric (cm)
    ID_LSM = 20/10; %cm
    ID_YSZ = ID_LSM + LSM_T;
    ID_NIYSZ = ID_YSZ + YSZ_T; 
    OD_NIYSZ = ID_NIYSZ + NIYSZ_T;
    LENGTH_TUBE = 1800/10;

    N_TUBES = 500;

    VOLUME_TANK    = 200*200*200; %cm3
    VOLUME_OD_TUBE = ((OD_NIYSZ/2)^2)*3.14159*LENGTH_TUBE + 2*3.14159*((OD_NIYSZ/2)^2);
    TOTAL_OD_VOL_TUBE = N_TUBES*VOLUME_OD_TUBE;

    VOLUME_ID_TUBE = ((ID_LSM/2)^2)*3.14159*LENGTH_TUBE + 2*3.14159*((ID_LSM/2)^2);
    TOTAL_ID_VOL_TUBE = N_TUBES*VOLUME_ID_TUBE;

    VOLUME_CATHODIC = VOLUME_TANK - TOTAL_OD_VOL_TUBE;
    VOLUME_ANODIC   = TOTAL_ID_VOL_TUBE;

    AREA_OD_TUBE = OD_NIYSZ*3.14159*LENGTH_TUBE + 2*3.14159*((OD_NIYSZ/2)^2);
    TOTAL_OD_AREA_TUBE = N_TUBES*AREA_OD_TUBE;

    %molar fractions of constituents 
    OXYGEN_CATHODE = 0.0; %assumption is that oxygen will be consumed immediately
    %after formation therefore oxygen percentage would be negligible 

    %chemical constants
    M_water = 18.01528; %g/mol
    M_hydro = 1.00794*2; %g/mol
    
    %temp_indiv  %K
    %press_indiv  %MPa
    
    B_Dissociation_indiv = zeros(999,1);
    B_Ionic_indiv        = zeros(999,1);
    KINETIC_P_indiv      = zeros(999,1);
    P_O2c_indiv          = zeros(999,1);
    J_indiv              = zeros(999,1);
    NONIDEAL_P_indiv     = zeros(999,1);
    IDEAL_P_indiv        = zeros(999,1);
    nohmic_a             = zeros(999,1);
    nact_a               = zeros(999,1);
    nconH2_a             = zeros(999,1);
    nconO2_a             = zeros(999,1);
    
    
    %Setting Variables
    %===============================================================================
    row = round((temp_indiv - 20 - 273.15)/5);
    col = round(press_indiv - 22);
    P_indiv              = POTENTIAL(col, row);
    %Calculating K Factor
    %===============================================================================
    K_indiv = exp(-GIBBS_FORMATION(col,row)*1000/(8.314*temp_indiv));
    
    %Electronic Region Check
    %===============================================================================
    Min_950 = 10^(-16);
    Max_950 = 10^7;
    Min_650 = 10^(-31);
    Max_650 = 10^12;
      
    Min_l = Min_650 + ((Min_950 - Min_650)/(950-650))*(temp_indiv - (650 + 273.15));
    Max_l = Max_650 + ((Max_950 - Max_650)/(950-650))*(temp_indiv - (650 + 273.15));
    
    E_Dissoc_indiv = -2.84215 + 5.02487*10^(-4)*temp_indiv -1.23698*10^(-8)*temp_indiv^2;
    
    for rxn_index = 1:1:999
      rxn = rxn_index/1000;
      
      %Mapping Thermodynamic Potential
      %===============================================================================
      [pot_ideal, non_ideal] = thermo__potentialb(rxn,OXYGEN_ANODE,temp_indiv,press_indiv,K_indiv);
      IDEAL_P_indiv   (rxn_index) = pot_ideal;
      NONIDEAL_P_indiv(rxn_index) = non_ideal;
      
      %Electronic Region Check
      %===============================================================================
      P_O2c_indiv(rxn_index) = sqrt(K_indiv)*(1 - rxn)/rxn;
      
      if(P_O2c_indiv(rxn_index) < Min_l || P_O2c_indiv(rxn_index) > Max_l)
          B_Ionic_indiv (rxn_index) = 0;
      else
          B_Ionic_indiv (rxn_index) = 1;
      end
      
      %Zirconia Dissociation Check
      %===============================================================================
      
      if(E_Dissoc_indiv < NONIDEAL_P_indiv(rxn_index))
          B_Dissociation_indiv (rxn_index) = 1;
      else
          B_Dissociation_indiv (rxn_index) = 0;
      end
      
      J_indiv(rxn_index) = 2.50; %A/cm2 current density intitial convergence value
      Pre_kinetic = NONIDEAL_P_indiv(rxn_index);
      KINETIC_P_indiv(rxn_index) = -100; %just putting a case that would initiate the while loop. 
      SF_Dissociation  = .9;              
      tol_C = .01; %in amps
      tol_V = .1;
      Top_J = 5.00;
      Bot_J = 0.00;
      criteria_kinetic_C = Top_J - Bot_J;
      criteria_kinetic_V = abs((abs(Pre_kinetic + KINETIC_P_indiv(rxn_index)) - abs(E_Dissoc_indiv*(SF_Dissociation))));
      while ((criteria_kinetic_C > tol_C) || (criteria_kinetic_V > tol_V))        
           %Mapping Kinetic Potential
           %===============================================================================
           [nohmic, nact, nconH2, nconO2, kinetic_pot] = kinetics(-J_indiv(rxn_index), ...
                               temp_indiv, press_indiv, rxn, ...
                               YSZ_T*10000, NIYSZ_T, LSM_T, ...
                               NIYSZ_Porosity, NIYSZ_Torosity, NIYSZ_radius, ...
                               LSM_Porosity, LSM_Torosity, LSM_radius, OXYGEN_ANODE);
            if(isnan(kinetic_pot))
                kinetic_pot = -100;
            end
            KINETIC_P_indiv(rxn_index) = kinetic_pot;
            nohmic_a(rxn_index)  = nohmic;
            nact_a(rxn_index)    = nact;
            nconH2_a(rxn_index)  = nconH2;
            nconO2_a(rxn_index)  = nconO2;
            
            criteria_kinetic_C = (Top_J - Bot_J);
            criteria_kinetic_V = abs((abs(Pre_kinetic + KINETIC_P_indiv(rxn_index)) - abs(E_Dissoc_indiv*(SF_Dissociation))));
            
            if((criteria_kinetic_C > tol_C) || (criteria_kinetic_V > tol_V)) 
                if(abs(Pre_kinetic + KINETIC_P_indiv(rxn_index)) > abs(E_Dissoc_indiv*(SF_Dissociation)))
                    Top_J = J_indiv(rxn_index);
                else 
                    Bot_J = J_indiv(rxn_index);
                end
            else
                 break
            end
            if ( J_indiv(rxn_index) > 4.999)
                 break;
            end
            J_indiv(rxn_index) = (Bot_J + Top_J)/2;
            
      end   
    end
%% Illustration of current + voltage subcomponents 
figure(20)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

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
plot(X_rxn, NONIDEAL_P_indiv + nohmic_a, '--r');
hold on
plot(X_rxn, NONIDEAL_P_indiv + nohmic_a + nact_a, '--b');
hold on
plot(X_rxn, NONIDEAL_P_indiv + nohmic_a + nact_a + nconH2_a, '--m');
hold on
plot(X_rxn, NONIDEAL_P_indiv + nohmic_a + nact_a + nconH2_a + nconO2_a, '--c');
hold on
l2 = plot(X_rxn,Zir, ':r');
set(l2,'linewidth',2);
hold on

xlabel('Reaction Completion [%]') % x-axis label
ylabel('Potential [V]') % y-axis label
set(gca,'XLim',[0 100])
set(gca,'XTick',(0:10:100))
set(gca,'YLim',[-3.5 0])
set(gca,'YTick',((-3.5):.5:0))
yyaxis right
l1 = plot(X_rxn,J_indiv,'-b');
set(l1,'linewidth',2);
ylabel('Current [A/cm^2]') % y-axis label
set(gca,'YLim',[0 max(J_indiv)*1.5])
set(gca,'YTick',(0:.5:max(J_indiv)*1.5))
legend('Thermodynamic','Non-Ideal Thermo Addition', 'Kinetic Addition',...
'Ohmic Resistance', 'Activation Polarization', 'Hydrogen Concentration Polar.',...
'Oxygen Concentration Polar.','Zirconia Dissociation','Current Density')
hold off

saveas(gcf,'Current&Voltage_Rxn_' + filename,'fig')
saveas(gcf,'Current&Voltage_Rxn_' + filename,'png')
%% Kinetic BreakDown
figure(21)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

C = horzcat(-nohmic_a,  -nact_a, -nconH2_a, -nconO2_a);
baseline = 0;
area(X_rxn,C,baseline);
legend('Ohmic Resistance', 'Activation Polarization', 'Hydrogen Concentration Polar.',...
'Oxygen Concentration Polar.')
xlabel('Reaction Completion [%]') % x-axis label
ylabel('Potential [V]') % y-axis label
set(gca,'XLim',[0 100])
set(gca,'XTick',(0:10:100))
saveas(gcf,'Kinetic_Voltage_Rxn_' + filename,'fig')
saveas(gcf,'Kinetic_Voltage_Rxn_' + filename,'png')

%==========================Beginning of Reaction

figure(22)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

area(X_rxn,C,baseline);
legend('Ohmic Resistance', 'Activation Polarization', 'Hydrogen Concentration Polar.',...
'Oxygen Concentration Polar.')
xlabel('Reaction Completion [%]') % x-axis label
ylabel('Potential [V]') % y-axis label
set(gca,'XLim',[0 1])
set(gca,'XTick',(0:.1:1))
saveas(gcf,'Kinetic_Voltage_RxnS_' + filename,'fig')
saveas(gcf,'Kinetic_Voltage_RxnS_' + filename,'png')

%==========================End of Reaction

figure(23)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

area(X_rxn,C,baseline);
legend('Ohmic Resistance', 'Activation Polarization', 'Hydrogen Concentration Polar.',...
'Oxygen Concentration Polar.')
xlabel('Reaction Completion [%]') % x-axis label
ylabel('Potential [V]') % y-axis label
set(gca,'XLim',[99 100])
set(gca,'XTick',(99:.1:100))
saveas(gcf,'Kinetic_Voltage_RxnE_' + filename,'fig')
saveas(gcf,'Kinetic_Voltage_RxnE_' + filename,'png')


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
           
           J2   = (-49+(J_index-1));
           rxn2 = ((rxn_index2-1)*49 + 1)/100;
           J_indiv2(J_index, rxn_index2) = J2;
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
            nact_a2  (J_index, rxn_index2)  = nact;
            nconH2_a2(J_index, rxn_index2)  = nconH2;
            nconO2_a2(J_index, rxn_index2)  = nconO2; 
            
    end 
end

%%
figure (997)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
s = 2; % 1=1%, 2=50%, 3=99%
Pot2       = zeros(100,1);
Zir2       = zeros(100,1);
for i = 1:length(Pot2) 
    % set each element to 0
    Pot2(i) = P_indiv;
    Zir2(i) = E_Dissoc_indiv;
end

title(['Voltage vs Current at ' , num2str(temp_indiv), ' K & ' , num2str(press_indiv) , ' MPa'])
plot(Pot2, J_indiv2(:,s),'-g');
hold on
plot(NONIDEAL_P_indiv2(:,s), J_indiv2(:,s), '-r');
hold on
plot(NONIDEAL_P_indiv2(:,s) + KINETIC_P_indiv2(:,s), J_indiv2(:,s), '-b');
hold on
plot(NONIDEAL_P_indiv2(:,s) + nohmic_a2(:,s), J_indiv2(:,s),  '--r');
hold on
plot(NONIDEAL_P_indiv2(:,s) + nohmic_a2(:,s) + nact_a2(:,s), J_indiv2(:,s),'--b');
hold on
plot(NONIDEAL_P_indiv2(:,s) + nohmic_a2(:,s) + nact_a2(:,s) + nconH2_a2(:,s), J_indiv2(:,s), '--m');
hold on
plot(NONIDEAL_P_indiv2(:,s) + nohmic_a2(:,s) + nact_a2(:,s) + nconH2_a2(:,s) + nconO2_a2(:,s), J_indiv2(:,s),'--c');
hold on
l2 = plot(Zir2, J_indiv2(:,s), ':r');
set(l2,'linewidth',2);
set(gca,'XLim',[-10 1])
set(gca,'XTick',((-10):1:1))
set(gca,'YLim',[-10 1])
set(gca,'YTick',((-10):1:1))
hold on
line(xlim(), [0,0], 'LineWidth', .5, 'Color', 'k');
hold on
xlabel('Potential [V]') % x-axis label
ylabel('Current [A/cm^2]') % y-axis label

legend({'Thermodynamic','Non-Ideal Thermo Addition', 'Kinetic Addition',...
'Ohmic Resistance', 'Activation Polarization', 'Hydrogen Concentration Polar.',...
'Oxygen Concentration Polar.','Zirconia Dissociation'},'Location','northwest')
hold off

saveas(gcf,'VoltagevCurrent_Whole_' + filename,'fig')
saveas(gcf,'VoltagevCurrent_Whole_' + filename,'png')
%%
figure (998)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
index_r = 1;
Pot2       = zeros(100,1);
Zir2       = zeros(100,1);
for i = 1:length(Pot2) 
    % set each element to 0
    Pot2(i) = P_indiv;
    Zir2(i) = E_Dissoc_indiv;
end

title(['Voltage vs Current at ' , num2str(temp_indiv), ' K & ' , num2str(press_indiv) , ' MPa'])

plot(nohmic_a2(:,index_r),       J_indiv2(:,index_r), '--r');
hold on
plot(nact_a2(:,index_r),         J_indiv2(:,index_r), '--b');
hold on
plot(nconH2_a2(:,index_r),       J_indiv2(:,index_r), '--m');
hold on
plot(nconO2_a2(:,index_r),       J_indiv2(:,index_r), '--c');
hold on
set(gca,'XLim',[-5 1])
set(gca,'XTick',((-5):1:1))
set(gca,'YLim',[-5 1])
set(gca,'YTick',((-5):1:1))
hold on
line(ylim(), [0,0], 'LineWidth', .5, 'Color', 'k');
%r = rectangle('Position',[-5,0,10,5],'FaceColor',[0 .5 .5]);
%facealpha(r,.1);
%hold on
xlabel('Potential [V]') % x-axis label
ylabel('Current [A/cm^2]') % y-axis label

legend({'Ohmic Resistance', 'Activation Polarization', 'Hydrogen Concentration Polar.',...
'Oxygen Concentration Polar.'},'Location','northwest')
hold off
saveas(gcf,'VoltagevCurrent_Comp_' + filename,'fig')
saveas(gcf,'VoltagevCurrent_Comp_' + filename,'png')
%% Reaction Completion Dependency
figure (999)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
Pot2       = zeros(100,1);
Zir2       = zeros(100,1);
for i = 1:length(Pot2) 
    % set each element to 0
    Pot2(i) = P_indiv;
    Zir2(i) = E_Dissoc_indiv;
end

title(['Voltage vs Current at ' , num2str(temp_indiv), ' K & ' , num2str(press_indiv) , ' MPa'])

plot(NONIDEAL_P_indiv2(:,1) + nohmic_a2(:,1) + nact_a2(:,1) + nconH2_a2(:,1) + nconO2_a2(:,1), J_indiv2(:,1),'--r');
hold on
plot(NONIDEAL_P_indiv2(:,2) + nohmic_a2(:,2) + nact_a2(:,2) + nconH2_a2(:,2) + nconO2_a2(:,2), J_indiv2(:,2),'--b');
hold on
plot(NONIDEAL_P_indiv2(:,3) + nohmic_a2(:,3) + nact_a2(:,3) + nconH2_a2(:,3) + nconO2_a2(:,3), J_indiv2(:,3),'--m');
hold on
l2 = plot(Zir2, J_indiv2(:,2), ':r');
set(l2,'linewidth',2);
set(gca,'XLim',[-5 1])
set(gca,'XTick',((-5):1:1))
set(gca,'YLim',[-5 1])
set(gca,'YTick',((-5):1:1))
hold on
line(xlim(), [0,0], 'LineWidth', .5, 'Color', 'k');
hold on
line(ylim(), [0,0], 'LineWidth', .5, 'Color', 'k');
hold on
xlabel('Potential [V]') % x-axis label
ylabel('Current [A/cm^2]') % y-axis label

legend('0.1% Reaction', '50.0% Reaction', '99.9% Reaction','Zirconia Dissociation')
hold off

saveas(gcf,'VoltagevCurrent_RxnV_' + filename,'fig')
saveas(gcf,'VoltagevCurrent_RxnV_' + filename,'png')
%% Energy Post Processing
%Start of Reaction 
rxn_initial = .05;
rxn_final   = .999;
%Calculation of Energy Produced by Reaction 
Final_Voltage = Pot + NONIDEAL_P_indiv + KINETIC_P_indiv;
Final_Current = J_indiv*TOTAL_OD_AREA_TUBE;

Final_Power = -Final_Current.*Final_Voltage;

figure(25)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
title(['Power at ' , num2str(temp_indiv), ' K & ' , num2str(press_indiv) , ' MPa'])
plot(X_rxn,Final_Power, '-r');
xlabel('Reaction Completion [%]') % x-axis label
ylabel('Power [J]') % y-axis label

hold on

density_indiv_w = W_DENSITY(col, row);
density_indiv_h = H_DENSITY(col, row); %mol/L

i_moles_water = (VOLUME_CATHODIC*(rxn_final-rxn_initial)/1000)*density_indiv_w; 
%starting moles of water at the beginning of the reaction

Reduced_Scope_Power = Final_Power((rxn_initial*1000):(rxn_final*1000));
Reduced_Scope_Rxn = X_rxn((rxn_initial*1000):(rxn_final*1000));
Reduced_Scope_Current = Final_Current((rxn_initial*1000):(rxn_final*1000));

Averaged_Power = trapz(Reduced_Scope_Rxn,Reduced_Scope_Power)/(rxn_final-rxn_initial);
Averaged_Current = trapz(Reduced_Scope_Rxn, Reduced_Scope_Current)/(rxn_final-rxn_initial);

%1 Ampere = 1 Coloumb/second
charge_e = 1.60217662*10^(-19); %coloumbs/electron

Time_Compl = (i_moles_water*6.022*10^23*2)/(Averaged_Current/charge_e);%seconds
Total_Energy = Averaged_Power*Time_Compl; %Joule
Normalized_Energy = Total_Energy/i_moles_water;

resultant_bottling_pressure = bottling_benefit(temp_indiv,press_indiv);
%MPa

txt = {'Time to Completion: ' + string(round(Time_Compl,3,'significant')) + ' secs',...
       'Total Energy: ' + string(round(Total_Energy,3,'significant'))+ ' Joule',...
       'Noramalized Energy: ' + string(round(Normalized_Energy,3,'significant')) + ' Joule/mol'};
text(10,.4*max(Final_Power),txt,'FontSize',12)
hold off
saveas(gcf,'Power_Rxn_' + filename,'fig')
saveas(gcf,'Power_Rxn_' + filename,'png')
%% Stress Analysis press_indiv

%Stress_State(temp_indiv, 50 , ID_LSM, ID_YSZ, ID_NIYSZ, OD_NIYSZ,0);

%% Calculate How Much Energy Required From Bringing Hydrogen Gas to Required Pressure
%====================================================================================
%Energy_Compression = -8.314*(273.15 + 25)*ln(0.101325/resultant_bottling_pressure);
%Results are in Joules
%This is the ideal amount of energy required to pressure hydrogen gas to
%its bottled form. If the conclusions to 

%% Calculate Heat Transfer Loss from Tank (basic heat transfer model)
%====================================================================================
%this will assume that the tank is a sphere assuming basic effeciencies
%expected from that size

close all; 
clear frame_h
save(filename)
end


