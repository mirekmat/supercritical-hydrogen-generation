function [] = Supercritical(LSM_T, YSZ_T, NIYSZ_T, OXYGEN_ANODE)
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


%% Setting Up Graphing 
x = 25:5:1000;
y = 23:1:400;

[X,Y] = meshgrid(x,y);
%% Thermodynamic Pass using NIST data
load('C:\Users\Matthew\Documents\Supercritical Water\Thermodynamics\thermodynamic_data.mat')
GIBBS_FORMATION = -2*W_GIBBS + (O_GIBBS + 2*H_GIBBS)+50.30007462*2 +237.141*2;
POTENTIAL = -GIBBS_FORMATION/(4*96.485);

%Thermodynamic Potential without considering effects of concentrations 
figure(1)
contour(X,Y,POTENTIAL,'ShowText','on')

xlabel('Temperature (°C)') % x-axis label
ylabel('Pressure (MPa)') % y-axis label

%% Clear Variables

%allocating space beforehand
B_Ionic              = zeros(378,196, 3);
NONIDEAL_P           = zeros(378,196);
IDEAL_P              = zeros(378,196);
P_O2c                = zeros(378,196, 3);
Bottle               = zeros(378,196);
Failure              = zeros(378,196);

%% Running Main Supercritical Scripts 

for col = 1:1:126
  for row = 1:1:196
            %Setting Variables
      %===============================================================================
      temp = (5*row + 20 + 273.15);
      press = (col + 22);
      
      %% Calculating K Factor
      %===============================================================================
      K_indiv = exp(-GIBBS_FORMATION(col,row)*1000/(8.314*temp));
     
      %% Stress Contour
      %===============================================================================
      Failure(col,row) = Stress_State_Sim(temp, press , ID_LSM, ID_YSZ, ID_NIYSZ, OD_NIYSZ,0); 
      %% Bottle Benefit
      %===============================================================================
      %Hydrogen Storage Tank 

      %Assuming that the tank is much larger 
      %Type I
      %Metal tank (steel/aluminum)
      %Approximate maximum pressure, aluminum 175 bars (17.5 MPa; 2,540 psi), steel 200 bars (20 MPa; 2,900 psi).
      Type1 = 20; %20 MPa
      %Type II
      %Metal tank (aluminum) with filament windings like glass fiber/aramid or carbon fiber around the metal cylinder.[6] See composite overwrapped pressure vessel.Approximate maximum pressure, aluminum/glass 263 bars (26.3 MPa; 3,810 psi), steel/carbon or aramide 299 bars (29.9 MPa; 4,340 psi).
      Type2 = 29.9; %20 MPa
      %Type III
      %Tanks made from composite material, fiberglass/aramid or carbon fiber with a metal liner (aluminum or steel). See metal matrix composite.
      %Approximate maximum pressure, aluminum/glass 305 bars (30.5 MPa; 4,420 psi), aluminum/aramide 438 bars (43.8 MPa; 6,350 psi), aluminium/ carbon 700 bars (70 MPa; 10,000 psi).
      Type3 = 70; %70 MPa
      
      Bottle(col,row) = bottling_benefit(temp, press);
          for rxn_index = 1:1:3
               rxn = ((rxn_index-1)*49 + 1)/100;
                   %% Electronic Region Check
                   %===============================================================================
                   Min_950 = 10^(-16);
                   Max_950 = 10^7;
                   Min_650 = 10^(-31);
                   Max_650 = 10^12;
      
                   Min_l = (Min_650 + ((Min_950 - Min_650)/(950-650))*(temp - (650 + 273.15)))/10;
                   Max_l = (Max_650 + ((Max_950 - Max_650)/(950-650))*(temp - (650 + 273.15)))/10;
                   t = (sqrt(K_indiv*(1 - rxn)/rxn))*press;
                   t2 = sqrt(K_indiv)*((1 - rxn)/rxn)*press;
                   P_O2c(col,row,rxn_index) = t;
      
                   if(P_O2c(col,row, rxn_index) > Min_l && P_O2c(col,row, rxn_index) < Max_l)
                       B_Ionic (col,row, rxn_index) = 1;
                   else
                       B_Ionic (col,row, rxn_index) = 0;
                   end 
                   
                  %% Mapping Thermodynamic Potential
                  %===============================================================================
                  [pot_ideal, non_ideal] = thermo__potentialb(rxn,OXYGEN_ANODE,temp,press,K_indiv);
                  IDEAL_P   (col,row, rxn_index) = pot_ideal;
                  NONIDEAL_P(col,row, rxn_index) = non_ideal;
               
          end
  end
end

save("PT_main")

%% Post Processing Data 
      %Partial Pressure
figure(20)

      contour(X,Y,B_Ionic(:,:,3),1,'ShowText','on')
      hold on
      contour(X,Y,B_Ionic(:,:,2),1,'ShowText','on')
      hold on
      contour(X,Y,B_Ionic(:,:,1),1,'ShowText','on')
      xlabel('Temperature (°C)') % x-axis label
      ylabel('Pressure (MPa)') % y-axis label
        
      hold off
      
saveas(gcf,'Ionic_Electronic','fig')
saveas(gcf,'Ionic_Electronic','png')
%% Failure of Membrane
      figure(25)
      contour(X,Y,Failure,9,'ShowText','on')
      hold on
      v = [.01,.01];
      contour(X,Y,Failure,v,'ShowText','on')
      hold on
      v = [.001,.001];
      contour(X,Y,Failure,v,'ShowText','on')
      xlabel('Temperature (°C)') % x-axis label
      ylabel('Pressure (MPa)') % y-axis label
      hold off
      
saveas(gcf,'Failure','fig')
saveas(gcf,'Failure','png')
%% Graphing Benefit
      figure(26)
      contour(X,Y,Bottle,[50,50],'ShowText','on')
      hold on
      contour(X,Y,Bottle,[100,100],'ShowText','on')
      hold on
      contour(X,Y,Bottle,[200,200],'ShowText','on')
      hold on
      v = [Type1,Type1];
      contour(X,Y,Bottle,v,'ShowText','on')
      hold on
      v = [Type2,Type2];
      contour(X,Y,Bottle,v,'ShowText','on')
      hold on
      v = [Type3,Type3];
      contour(X,Y,Bottle,v,'ShowText','on')
      title('Bottling Benefit')
      xlabel('Temperature (°C)') % x-axis label
      ylabel('Pressure (MPa)') % y-axis label
      hold off
      
saveas(gcf,'Bottling_Figure','fig')
saveas(gcf,'Bottling_Figure','png')

%% Graphing Benefit vs Failure 
      figure(99)
      contour(X,Y,Failure,9,'ShowText','on')
      hold on
      v = [.01,.01];
      contour(X,Y,Failure,v,'ShowText','on')
      hold on
      v = [.001,.001];
      contour(X,Y,Failure,v,'ShowText','on')
      hold on
      v = [Type1,Type1];
      contour(X,Y,Bottle,v,'ShowText','on')
      hold on
      v = [Type2,Type2];
      contour(X,Y,Bottle,v,'ShowText','on')
      hold on
      v = [Type3,Type3];
      contour(X,Y,Bottle,v,'ShowText','on')
      title('Bottling Benefit')
      xlabel('Temperature (°C)') % x-axis label
      ylabel('Pressure (MPa)') % y-axis label
      hold off
      
saveas(gcf,'FailureBottle_Figure','fig')
saveas(gcf,'FailureBottle_Figure','png')
end






      










