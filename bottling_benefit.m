function [Final] = bottling_benefit(To,Po)

%% Hydrogen Storage Tank 

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


%H2 Fugacity Data & Calculations
H2_Tc = 33.19;        %K
H2_Pc = 1.31;         %MPa
H2_Vc = 6.42*10^-5;   %m3/mol
H2_w  = -.2320;       %w

alpha_H2  = (1+(0.37464 + 1.5422*H2_w - 0.26992*H2_w^2)*(1-sqrt(To/H2_Tc)))^2;
A_H2  = 0.45726*(8.314^2)*(H2_Tc)^2*alpha_H2/(H2_Pc*10^6);
B_H2  = 0.0778*8.3154*H2_Tc/(H2_Pc*10^6);
%-------------------------------------------------
%Calculation at variable pressure and temperature
A1 =  A_H2*Po*10^6/(8.314^2*To^2);
B1 =  B_H2*Po*10^6/(8.314*To);

p1 = [1 -(1-B1) (A1-2*B1-3*B1^2) -(A1*B1-B1^2-B1^3)];

r1 = roots(p1);
Z1 = max(r1);

%Calculations at 273.15 K & 30 MPa
Ts = 273.15 + 25;
Ps = 50; % initial guess
P  = 1;  % initialization

while (abs(Ps - P) > .01) 
    
    A2 =  A_H2*Ps*10^6/(8.314^2*Ts^2);
    B2 =  B_H2*Ps*10^6/(8.314*Ts);

    p2 = [1 -(1-B2) (A2-2*B2-3*B2^2) -(A2*B2-B2^2-B2^3)];

    r2 = roots(p2);
    Z2 = max(r2);

    P  = Po*(Ts*Z2)/(To*Z1);
    
    Ps = Ps - (Ps - P)/2;
    
end

Final = P;

%Output is the final pressure after cooling down the supercritical solution
%to room temperature. This is under the assumption that no more work needs
%to be added to the system. Therefore a basic baseline curve is created of
%zones that cannot be used in the supercrtical region of hydrogen gas since
%it violates the condition that the compression of gas would soley be
%derived from the thermal energy that the supercritical region originally
%holds. 

