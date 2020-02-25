function [pot_ideal, non_ideal] = thermo__potentialb(H2_y,O2_ya, To, Po, K)
%Calls for fugacity data for both the anodic and cathodic side and
%calculates the concentration potentional. Also displays potential for
%ideal case. Anode side is only oxygen and argon gas as the inert
%material. This not a function of reaction progress but is regulated. 

%Calculating volume fraction of 
Ar_y = 1 - O2_ya;
O2_yc = 0;
H2O_y = 1 - H2_y;


[H2O_theta,H2_theta,O2_theta_c, Zc] = fugacity_cathode(H2O_y,H2_y,O2_yc,To, Po);
[Ar_theta,O2_theta_a, Za]           = fugacity_anode(Ar_y,O2_ya,To, Po);

f_w  = H2O_theta*H2O_y*Po*10^6;
f_h  = H2_theta*H2_y*Po*10^6;
f_oa = O2_theta_a*O2_ya*Po*10^6;

ID = (K*(H2O_y*Po*10^6)^2)/(O2_ya*Po*10^6*(H2_y*Po*10^6)^2);
LNID = reallog(ID);
pot_ideal = (8.314*To*LNID/(96485*4));

NID = (K*f_w^2)/(f_oa*f_h^2);
LNNID = reallog(NID);
non_ideal = (8.314*To*LNNID/(96485*4));


end