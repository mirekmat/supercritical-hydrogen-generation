function [O2_D] = diffusion_oxygen(To,Po)

Dp = 9.77*10^(-6)*To^1.736;
%constants use pressure in atm 
%hydrogen diffusion is supercritical water requires a mixed density while
%oxygen in argon requires pressure of the system

%1 MPa = 9.869 atm
Patm = Po*9.869;
O2_D  = Dp/Patm;

%cm^2/s diffusion unit
end