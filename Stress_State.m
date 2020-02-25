function [PROB_FAILURE] = Stress_State(temp_indiv, press_indiv, ID_LSM, ID_YSZ, ID_NIYSZ, OD_NIYSZ, PERTUR)

%were added on 
TEMP = temp_indiv;
%replace later with the actual pressure amount of the system, they should
%be equal, but useful in understanding what pressure perturbations do to a
%system. 
P_inner = press_indiv*1000000 + PERTUR*1000000/2;
P_outer = press_indiv*1000000 - PERTUR*1000000/2;

PVD_temp = temp_indiv;  %273.15 + 800; %temperature that layers of YSZ, Ni- NYSZ, LSM 

%%Import Data from Tables for Material Constants
%Original data is from a

%%
load('material.mat')

CTE_LSM_x = table2array(LSM_CTE(:,1));
CTE_YSZ_x = table2array(YSZ_CTE(:,1));
CTE_NIY_x = table2array(NIYSZ_CTE(:,1));
E_LSM_x = table2array(LSM_E(:,1));
E_YSZ_x = table2array(YSZ_E(:,1));
E_NIY_x = table2array(NIYSZ_E(:,1));
CTE_LSM_y = table2array(LSM_CTE(:,2));
CTE_YSZ_y = table2array(YSZ_CTE(:,2));
CTE_NIY_y = table2array(NIYSZ_CTE(:,2));
E_LSM_y = table2array(LSM_E(:,2));
E_YSZ_y = table2array(YSZ_E(:,2));
E_NIY_y = table2array(NIYSZ_E(:,2));

%Eliminates any spaces from the arrays, since data sets are uneven.
CTE_LSM_x( all( ~any( CTE_LSM_x), 1 ), : ) = []; 
CTE_YSZ_x( all( ~any( CTE_YSZ_x), 1 ), : ) = []; 
CTE_NIY_x( all( ~any( CTE_NIY_x), 1 ), : ) = []; 
E_LSM_x( all( ~any( E_LSM_x), 1 ), : )     = []; 
E_YSZ_x( all( ~any( E_YSZ_x), 1 ), : )     = []; 
E_NIY_x( all( ~any( E_NIY_x), 1 ), : )     = [];
CTE_LSM_y( all( ~any( CTE_LSM_y), 1 ), : ) = []; 
CTE_YSZ_y( all( ~any( CTE_YSZ_y), 1 ), : ) = []; 
CTE_NIY_y( all( ~any( CTE_NIY_y), 1 ), : ) = []; 
E_LSM_y( all( ~any( E_LSM_y), 1 ), : )     = []; 
E_YSZ_y( all( ~any( E_YSZ_y), 1 ), : )     = []; 
E_NIY_y( all( ~any( E_NIY_y), 1 ), : )     = []; 

%
CTE_LSM_x(isnan(CTE_LSM_x))=[];
CTE_LSM_y(isnan(CTE_LSM_y))=[];
CTE_YSZ_x(isnan(CTE_YSZ_x))=[];
CTE_YSZ_y(isnan(CTE_YSZ_y))=[];
CTE_NIY_x(isnan(CTE_NIY_x))=[];
CTE_NIY_y(isnan(CTE_NIY_y))=[];

E_LSM_x(isnan(E_LSM_x))=[];
E_LSM_y(isnan(E_LSM_y))=[];
E_YSZ_x(isnan(E_YSZ_x))=[];
E_YSZ_y(isnan(E_YSZ_y))=[];
E_NIY_x(isnan(E_NIY_x))=[];
E_NIY_y(isnan(E_NIY_y))=[];
%Constant extrapolation to cover both end points, spline curve
%extrapolation didn't seem to allow realistic values for larger
%perturbations of temperature. This seemed to the most conservative
%approximation. 
CTE_LSM_x = [200; CTE_LSM_x; 1500];
CTE_YSZ_x = [200; CTE_YSZ_x; 1500];
CTE_NIY_x = [200; CTE_NIY_x; 1500];
E_LSM_x = [200; E_LSM_x; 1500];
E_YSZ_x = [200; E_YSZ_x; 1500];
E_NIY_x = [200; E_NIY_x; 1500];

CTE_LSM_y = [CTE_LSM_y(1); CTE_LSM_y; CTE_LSM_y(length(CTE_LSM_y))];
CTE_YSZ_y = [CTE_YSZ_y(1); CTE_YSZ_y; CTE_YSZ_y(length(CTE_YSZ_y))];
CTE_NIY_y = [CTE_NIY_y(1); CTE_NIY_y; CTE_NIY_y(length(CTE_NIY_y))];

E_LSM_y = [E_LSM_y(1); E_LSM_y; E_LSM_y(length(E_LSM_y))];
E_YSZ_y = [E_YSZ_y(1); E_YSZ_y; E_YSZ_y(length(E_YSZ_y))];
E_NIY_y = [E_NIY_y(1); E_NIY_y; E_NIY_y(length(E_NIY_y))];


%%

%regression model ELASTIC MODULUS 

           
%ftE1 = fittype('a-b*T*exp(-c/T) + d*(T-e*c + abs(T -e*c))*exp(-c/T)',...
%             'options',foEY,'independent', 'T' ,'dependent', 'E');
%ftE = fittype('a+b*T+c*(T-g)^2+d*(T-h)^3+e*(T-f)^4+f*T^5',...
%             'options',foEY,'independent', 'T' ,'dependent', 'E');

figure(79)
         
plot(E_YSZ_x,E_YSZ_y,'o')
hold on
[curve_YSZ_E,gof_YSZ_E] = fit(E_YSZ_x,E_YSZ_y, 'smoothingspline','SmoothingParam',.001);
plot(curve_YSZ_E,'m')
hold on


plot(E_NIY_x,E_NIY_y,'o')
hold on
[curve_NIY_E,gof_NIY_E] = fit(E_NIY_x,E_NIY_y,'smoothingspline','SmoothingParam',.001);
plot(curve_NIY_E,'m')
hold on


plot(E_LSM_x,E_LSM_y,'o')
hold on
[curve_LSM_E,gof_LSM_E] = fit(E_LSM_x,E_LSM_y,'smoothingspline','SmoothingParam',.001);
plot(curve_LSM_E,'m')
hold on

legend('YSZ Data','YSZ Fit','NIY Data','NIY Fit','LSM Data','LSM Fit')
xlabel("Temperature [K]");
ylabel("Elastic Modulus [GPa]");

hold off

saveas(gcf,'E_Fitted Data' + filename,'fig')
saveas(gcf,'E_Fitted Data' + filename,'png')


figure(89)

plot(CTE_YSZ_x,CTE_YSZ_y,'o')
hold on
[curve_YSZ_CTE,gof_YSZ_CTE] = fit(CTE_YSZ_x,CTE_YSZ_y, 'smoothingspline','SmoothingParam',.001);
plot(curve_YSZ_CTE,'m')
hold on


plot(CTE_NIY_x,CTE_NIY_y,'o')
hold on
[curve_NIY_CTE,gof_NIY_CTE] = fit(CTE_NIY_x,CTE_NIY_y,'smoothingspline','SmoothingParam',.001);
plot(curve_NIY_CTE,'m')
hold on


plot(CTE_LSM_x,CTE_LSM_y,'o')
hold on
[curve_LSM_CTE,gof_LSM_CTE] = fit(CTE_LSM_x,CTE_LSM_y,'smoothingspline','SmoothingParam',.001);
plot(curve_LSM_CTE,'m')
hold on

legend('YSZ Data','YSZ Fit','NIY Data','NIY Fit','LSM Data','LSM Fit')
xlabel("Temperature [K]");
ylabel("CTE * 10^-6 [m/m]");
hold off


saveas(gcf,'CTE_Fitted Data' + filename,'fig')
saveas(gcf,'CTE_Fitted Data' + filename,'png')

%%

%calculating the CTE average over the temperature span that the memebrane
%is being exposed to, this will give a more accurate estimate of the
%thermal induced stressed to the system.
if(PVD_temp > TEMP)
     temp_scale  = TEMP:.01:PVD_temp;
elseif (PVD_temp < TEMP)
     temp_scale =  PVD_temp:.01:TEMP;
else 
     temp_scale = 0;
end 
int_LSM = integrate(curve_LSM_CTE,temp_scale,0)/(TEMP-PVD_temp);
growth_LSM = (max(int_LSM) - min(int_LSM))*10^-6;
int_YSZ = integrate(curve_YSZ_CTE,temp_scale,0)/(TEMP-PVD_temp);
growth_YSZ = (max(int_YSZ) - min(int_YSZ))*10^-6;
int_NIY = integrate(curve_NIY_CTE,temp_scale,0)/(TEMP-PVD_temp);
growth_NIY = (max(int_NIY) - min(int_NIY))*10^-6;
delta_T = (TEMP-PVD_temp);

E_NIY = 1000000000*curve_NIY_E(TEMP); %1000X to make from GPa to N/m2
E_LSM = 1000000000*curve_LSM_E(TEMP);
E_YSZ = 1000000000*curve_YSZ_E(TEMP);

%%
%Defining variables. 
B_1 = ID_LSM/100;
B_2 = ID_YSZ/100;
B_3 = ID_NIYSZ/100;
B_4 = OD_NIYSZ/100;

v_LSM = 0.3;
v_YSZ = 0.3;
v_NIY = 0.3;

%factor calculated for the elongation due to coeffecient of thermal
%expansion

if(PVD_temp == TEMP)
    growth_LSM = 0;
    growth_YSZ = 0;
    growth_NIY = 0;
end
    
A_LSM = 1/(1 - v_LSM)*(B_2^2 - B_1^2)*(growth_LSM);
A_YSZ = 1/(1 - v_YSZ)*(B_3^2 - B_2^2)*(growth_YSZ);
A_NIY = 1/(1 - v_NIY)*(B_4^2 - B_3^2)*(growth_NIY);

%%
syms C_LSM D_LSM C_YSZ D_YSZ C_NIY D_NIY

%Equilibrium Equation 1
eqn1 = (1/(1-2*v_LSM))*C_LSM - D_LSM/B_1^2 == -P_inner*(1 + v_LSM)/E_LSM;
%Equilibrium Equation 2
eqn2 = (1/(1-2*v_NIY))*C_NIY - D_NIY/B_4^2 == (1 + v_LSM)*(A_NIY/B_4^2 - P_outer/E_NIY);
%Equilibrium Equation 3
eqn3 = (1/(1-2*v_LSM))*C_LSM - D_LSM/B_2^2 - E_YSZ/(E_LSM*(1-2*v_LSM))*C_YSZ + E_YSZ*D_YSZ/(E_LSM*B_2^2) == (1+v_LSM)*A_LSM/B_2^2;
%Equilibrium Equation 4
eqn4 = (1/(1-2*v_YSZ))*C_YSZ - D_YSZ/B_3^2 - E_NIY/(E_YSZ*(1-2*v_YSZ))*C_NIY + E_NIY*D_NIY/(E_YSZ*B_3^2) == (1+v_YSZ)*A_YSZ/B_3^2;
%Equilibrium Equation 5
eqn5 = C_LSM + D_LSM/B_2^2 - C_YSZ - D_YSZ/B_2^2 == -(v_LSM)*(A_LSM/B_2^2);
%Equilibrium Equation 6
eqn6 = C_YSZ + D_YSZ/B_3^2 - C_NIY - D_NIY/B_3^2 == -(v_YSZ)*(A_YSZ/B_3^2);


%Solve the system of equations using solve. The inputs to solve are a 
%vector of equations, and a vector of variables to solve the equations for.
sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], ...
    [C_LSM, D_LSM, C_YSZ, D_YSZ, C_NIY, D_NIY]);

C_LSM_Sol = double(sol.C_LSM);
D_LSM_Sol = double(sol.D_LSM);
C_YSZ_Sol = double(sol.C_YSZ);
D_YSZ_Sol = double(sol.D_YSZ);
C_NIY_Sol = double(sol.C_NIY);
D_NIY_Sol = double(sol.D_NIY);

%computation of the stresses for each layer for the inner and outer portion
%of the layer

%LSM [Layer 1]
r_LSM = B_1:(B_2-B_1)/100:B_2;
sigma_r_L = (-E_LSM./((1 - v_LSM).*r_LSM.^2).*(r_LSM.^2 - B_1.^2).*(growth_LSM.*delta_T)...
          +(E_LSM./(1 + v_LSM)).*(C_LSM_Sol./(1-2*v_LSM) - D_LSM_Sol./r_LSM));
sigma_a_L = (-(growth_LSM*delta_T).*E_LSM/(1-v_LSM) + ...
          2*v_LSM*E_LSM*C_LSM_Sol./((1+v_LSM).*(1-2.*v_LSM)));
sigma_h_L = (E_LSM./((1-v_LSM).*r_LSM.^2).*(r_LSM.^2 - B_1.^2).*(growth_LSM.*delta_T)...
          -E_LSM.*(growth_LSM.*delta_T)./(1-v_LSM) + E_LSM.*(C_LSM_Sol./(1-2.*v_LSM)...
        - D_LSM_Sol./r_LSM)./(1+v_LSM));

%YSZ [Layer 2]
r_YSZ = B_2:(B_3-B_2)/100:B_3;
sigma_r_Y = (-E_YSZ./((1 - v_YSZ).*r_YSZ.^2).*(r_YSZ.^2 - B_2.^2).*(growth_YSZ.*delta_T)...
          +(E_YSZ./(1 + v_YSZ)).*(C_YSZ_Sol./(1-2.*v_YSZ) - D_YSZ_Sol./r_YSZ));
sigma_a_Y = (-(growth_YSZ.*delta_T).*E_YSZ./(1-v_YSZ) + ...
          2.*v_YSZ.*E_YSZ.*C_YSZ_Sol./((1+v_YSZ)*(1-2.*v_YSZ)));
sigma_h_Y = (E_YSZ./((1-v_YSZ).*r_YSZ.^2).*(r_YSZ.^2 - B_2.^2).*(growth_YSZ.*delta_T)...
          -E_YSZ.*(growth_YSZ.*delta_T)./(1-v_YSZ) + E_YSZ.*(C_YSZ_Sol./(1-2.*v_YSZ)...
        - D_YSZ_Sol./r_YSZ)./(1+v_YSZ));


%Ni-YSZ [Layer 3]
r_NIY = B_3:(B_4-B_3)/100:B_4;
sigma_r_N = (-E_NIY./((1 - v_NIY).*r_NIY.^2).*(r_NIY.^2 - B_3.^2).*(growth_NIY.*delta_T)...
          +(E_NIY./(1 + v_NIY)).*(C_NIY_Sol./(1-2.*v_NIY) - D_NIY_Sol./r_NIY));
sigma_a_N = (-(growth_NIY.*delta_T).*E_NIY/(1-v_NIY) + ...
          2.*v_NIY.*E_NIY.*C_NIY_Sol/((1+v_NIY)*(1-2.*v_NIY)));
sigma_h_N = (E_NIY./((1-v_NIY).*r_NIY.^2).*(r_NIY.^2 - B_3.^2)*(growth_NIY.*delta_T)...
          -E_NIY.*(growth_NIY.*delta_T)./(1-v_NIY) + E_NIY.*(C_NIY_Sol./(1-2.*v_NIY)...
        - D_LSM_Sol./r_NIY)./(1+v_NIY));
    
r_whole = cat(2,r_LSM, r_YSZ, r_NIY);

q = repelem(sigma_a_L,101);
w = repelem(sigma_a_Y,101);
e = repelem(sigma_a_N,101);

sigma_r_whole = cat(2,sigma_r_L, sigma_r_Y, sigma_r_N);
sigma_a_whole = cat(2,q,w,e);
sigma_h_whole = cat(2,sigma_h_L, sigma_h_Y, sigma_h_N);
top_limit = max([sigma_r_whole sigma_a_whole sigma_h_whole])/10^6;
bottom_limit = min([sigma_r_whole sigma_a_whole sigma_h_whole])/10^6;

figure (49)
plot (sigma_r_whole/10^6, '-r*');
hold on
plot (sigma_a_whole/10^6, '-b*');
hold on
plot (sigma_h_whole/10^6, '-go');
hold on 

line([101 101], [bottom_limit top_limit]);
hold on 
line([202 202], [bottom_limit top_limit]); 
hold off
legend('Radial Stress','Axial Stress','Hoop Stress')
xlabel("Radial Distance [normalized]");
ylabel("Stress [MPa]");

%% Weibulll Modulus
             %sigo [MPa] sigu [MPa]   m
LSM_298K  = [     52*1.50         0        6.7];% increase in best chance 
LSM_1073K = [     75*1.50         0        3.7];

NIY_298K  = [    115.0       0       13.2];
NIY_1073K = [     84.0       0       13.2];

YSZ_298K  = [    932.1       0       10.8];

LSM_Weibull = LSM_298K + (LSM_1073K-LSM_298K)*(temp_indiv - 273)/(1073-273);
NIY_Weibull = NIY_298K + (NIY_1073K-NIY_298K)*(temp_indiv - 273)/(1073-273);
YSZ_Weibull = YSZ_298K;

LSM_sigma_max = max([abs(sigma_r_L)  abs(sigma_a_L) abs(sigma_h_L)])/10^6;
NIY_sigma_max = max([abs(sigma_r_N)  abs(sigma_a_N) abs(sigma_h_N)])/10^6;
YSZ_sigma_max = max([abs(sigma_r_Y)  abs(sigma_a_Y) abs(sigma_h_Y)])/10^6;

LSM_F = 1-exp(-((LSM_sigma_max-LSM_Weibull(2))/LSM_Weibull(1))^LSM_Weibull(3));
NIY_F = 1-exp(-((NIY_sigma_max-NIY_Weibull(2))/NIY_Weibull(1))^NIY_Weibull(3));
YSZ_F = 1-exp(-((YSZ_sigma_max-YSZ_Weibull(2))/YSZ_Weibull(1))^YSZ_Weibull(3));

PROB_FAILURE = max([LSM_F NIY_F YSZ_F]);
end