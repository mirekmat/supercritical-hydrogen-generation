
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

%% Setting Up Graphing 
x = 25:5:1000;
y = 23:1:400;

[X,Y] = meshgrid(x,y);

for Tox = 25:5:1000
     for Pox = 23:1:400
         bottle(Pox-22, (Tox-20)/5) = bottling_benefit(Tox + 273.15, Pox); 
     end
end

%Thermodynamic Potential without considering effects of concentrations 
figure(91)
contour(X,Y,bottle,'ShowText','on')

xlabel('Temperature (°C)') % x-axis label
ylabel('Pressure (MPa)') % y-axis label

saveas(gcf,'Bottling_Benefit' + filename,'fig')
saveas(gcf,'Bottling_Benefit' + filename,'png')