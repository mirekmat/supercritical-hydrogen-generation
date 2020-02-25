% Pressure & Temperature Sensitivity Test Cases
                   

%Temp [K] ||Press [MPa] || LSM_T [cm] || YSZ_T [cm] ||  NIYSZ_T [cm]  || Anode O2 [-]
%Supercritical_indiv(273 + 700  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 700  ,    50       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 700  ,    100      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 700  ,    200      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 700  ,    400      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );

%Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 800  ,    50       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 800  ,    100      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 800  ,    200      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 800  ,    400      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );

%Supercritical_indiv(273 + 900  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 900  ,    50       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 900  ,    100      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 900  ,    200      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 900  ,    400      ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );

%% Material Thickness Sensitivity Test Cases
                    %Temp [K] ||Press [MPa] || LSM_T [cm] || YSZ_T [cm] ||  NIYSZ_T [cm]  || Anode O2 [-]
%Supercritical_indiv(273 + 800  ,    25       ,     1/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 800  ,    25       ,     3/10    ,   0.040/10  ,    0.150/10     ,    0.9      );

Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.030/10  ,    0.150/10     ,    0.9      );
%Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.050/10  ,    0.150/10     ,    0.9      );

Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.140/10     ,    0.9      );
%Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );
Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.160/10     ,    0.9      );

%% Oxygen Partial PressureSensitivity Test Cases
                    %Temp [K] ||Press [MPa] || LSM_T [cm] || YSZ_T [cm] ||  NIYSZ_T [cm]  || Anode O2 [-]
Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.1      );
%Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.3      );
Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.5      );
%Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.7      );
Supercritical_indiv(273 + 800  ,    25       ,     2/10    ,   0.040/10  ,    0.150/10     ,    0.9      );

%% Feasability of Supercritical Water
             % LSM_T [cm] || YSZ_T [cm]  ||  NIYSZ_T [cm]  || Anode O2 [-]
%Supercritical(  2/10       ,   0.040/10  ,    0.150/10     ,    0.9      );
