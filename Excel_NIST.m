%Script that utilizes excel files to obtain thermodynamic and fluid dynamic
%data for further calculations. Initializes variables that will be used by
%other functions. 

%% Import data from spreadsheet
%%

%% IMPORTING CONDENSED WATER DATA
WATER = zeros(404,14);
for i = 5:200
    
    j = i*5;
    %% Import the data
    [~, ~, raw] = xlsread('C:\Users\Matthew\Documents\Supercritical Water\Thermodynamics\Water.xlsx',int2str(j),'A2:N405');
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

    %% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells

    %% Create output variable
    Super = reshape([raw{:}],size(raw));

    %% Clear temporary variables
    clearvars raw R;
    
    WATER(:,:,i-4) = Super;
    
end

%%
%%IMPORTING HYDROGEN DATA
HYDROGEN = zeros(404,14);
for i = 5:200
    
    j = i*5;
    %% Import the data
    [~, ~, raw] = xlsread('C:\Users\Matthew\Documents\Supercritical Water\Thermodynamics\Hydrogen.xlsx',int2str(j),'A2:N405');
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

    %% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells

    %% Create output variable
    Super = reshape([raw{:}],size(raw));

    %% Clear temporary variables
    clearvars raw R;
    
    HYDROGEN(:,:,i-4) = Super;
    
end
%%

%%IMPORTING OXYGEN DATA
OXYGEN = zeros(404,14);
for i = 5:200
    
    j = i*5;
    %% Import the data
    [~, ~, raw] = xlsread('C:\Users\Matthew\Documents\Supercritical Water\Thermodynamics\Oxygen.xlsx',int2str(j),'A2:N405');
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

    %% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells

    %% Create output variable
    Super = reshape([raw{:}],size(raw));

    %% Clear temporary variables
    clearvars raw R;
    
    OXYGEN(:,:,i-4) = Super;
    
end
%%
%%Manipulations Water Graphing
W_TEMP = squeeze(WATER(:,1,:)) + 273.15; %Temperature (Kelvin )
W_DENSITY =  squeeze(WATER(:,3,:)); %Density (mol/l)	
W_VOLUME =  squeeze(WATER(:,4,:)); %Volume (l/mol)	
W_INTERNAL_ENERGY =   squeeze(WATER(:,5,:));%Internal Energy (kJ/mol)	
W_ENTHALPY =  squeeze(WATER(:,6,:)); %Enthalpy (kJ/mol)	
W_ENTROPY =  squeeze(WATER(:,7,:)); %Entropy (J/mol*K)	
W_CV =  squeeze(WATER(:,8,:)); %Cv (J/mol*K)	
W_CP =  squeeze(WATER(:,9,:)); %Cp (J/mol*K)	
W_SOUND_SPD =  squeeze(WATER(:,10,:)); %Sound Spd. (m/s)	
W_JOULE_THOMPSON =  squeeze(WATER(:,11,:));  %Joule-Thomson (K/MPa)	
W_VISCOSITY =  squeeze(WATER(:,12,:));  %Viscosity (Pa*s)	
W_THERMAL_CONDUCTIVITY =  squeeze(WATER(:,13,:)); %Therm. Cond. (W/m*K)	
W_PHASE =  squeeze(WATER(:,14,:));%Phase

%%
%%Isolating Supercrtical Region
W_TEMP(1:25,:)=[];
W_DENSITY(1:25,:) =  [];
W_VOLUME(1:25,:) =  [];	
W_INTERNAL_ENERGY(1:25,:) =  [];  
W_ENTHALPY(1:25,:) =  [];	
W_ENTROPY(1:25,:) =  [];
W_CV(1:25,:) =  [];
W_CP(1:25,:) =  [];
W_SOUND_SPD(1:25,:) =  []; 
W_JOULE_THOMPSON(1:25,:) =  [];
W_VISCOSITY(1:25,:) =  [];
W_THERMAL_CONDUCTIVITY(1:25,:) =  [];
W_PHASE(1:25,:) =  [];

W_TEMP(379,:)=[];
W_DENSITY(379,:) =  [];
W_VOLUME(379,:) =  [];	
W_INTERNAL_ENERGY(379,:) =  [];  
W_ENTHALPY(379,:) =  [];	
W_ENTROPY(379,:) =  [];
W_CV(379,:) =  [];
W_CP(379,:) =  [];
W_SOUND_SPD(379,:) =  []; 
W_JOULE_THOMPSON(379,:) =  [];
W_VISCOSITY(379,:) =  [];
W_THERMAL_CONDUCTIVITY(379,:) =  [];
W_PHASE(379,:) =  [];

%%
%%Manipulations Hydrogen Graphing

H_DENSITY =  squeeze(HYDROGEN(:,3,:)); %Density (mol/l)	
H_VOLUME =  squeeze(HYDROGEN(:,4,:)); %Volume (l/mol)	
H_INTERNAL_ENERGY =   squeeze(HYDROGEN(:,5,:));%Internal Energy (kJ/mol)	
H_ENTHALPY =  squeeze(HYDROGEN(:,6,:)); %Enthalpy (kJ/mol)	
H_ENTROPY =  squeeze(HYDROGEN(:,7,:)); %Entropy (J/mol*K)	
H_CV =  squeeze(HYDROGEN(:,8,:)); %Cv (J/mol*K)	
H_CP =  squeeze(HYDROGEN(:,9,:)); %Cp (J/mol*K)	
H_SOUND_SPD =  squeeze(HYDROGEN(:,10,:)); %Sound Spd. (m/s)	
H_JOULE_THOMPSON =  squeeze(HYDROGEN(:,11,:));  %Joule-Thomson (K/MPa)	
H_VISCOSITY =  squeeze(HYDROGEN(:,12,:));  %Viscosity (Pa*s)	
H_THERMAL_CONDUCTIVITY =  squeeze(HYDROGEN(:,13,:)); %Therm. Cond. (W/m*K)	
H_PHASE =  squeeze(HYDROGEN(:,14,:));%Phase

%%
%%Isolating Supercrtical Region
H_DENSITY(1:25,:) =  [];
H_VOLUME(1:25,:) =  [];	
H_INTERNAL_ENERGY(1:25,:) =  [];  
H_ENTHALPY(1:25,:) =  [];	
H_ENTROPY(1:25,:) =  [];
H_CV(1:25,:) =  [];
H_CP(1:25,:) =  [];
H_SOUND_SPD(1:25,:) =  []; 
H_JOULE_THOMPSON(1:25,:) =  [];
H_VISCOSITY(1:25,:) =  [];
H_THERMAL_CONDUCTIVITY(1:25,:) =  [];
H_PHASE(1:25,:) =  [];

H_DENSITY(379,:) =  [];
H_VOLUME(379,:) =  [];	
H_INTERNAL_ENERGY(379,:) =  [];  
H_ENTHALPY(379,:) =  [];	
H_ENTROPY(379,:) =  [];
H_CV(379,:) =  [];
H_CP(379,:) =  [];
H_SOUND_SPD(379,:) =  []; 
H_JOULE_THOMPSON(379,:) =  [];
H_VISCOSITY(379,:) =  [];
H_THERMAL_CONDUCTIVITY(379,:) =  [];
H_PHASE(379,:) =  [];

%%
%%Manipulations Oxygen Graphing

O_DENSITY =  squeeze(OXYGEN(:,3,:)); %Density (mol/l)	
O_VOLUME =  squeeze(OXYGEN(:,4,:)); %Volume (l/mol)	
O_INTERNAL_ENERGY =   squeeze(OXYGEN(:,5,:));%Internal Energy (kJ/mol)	
O_ENTHALPY =  squeeze(OXYGEN(:,6,:)); %Enthalpy (kJ/mol)	
O_ENTROPY =  squeeze(OXYGEN(:,7,:)); %Entropy (J/mol*K)	
O_CV =  squeeze(OXYGEN(:,8,:)); %Cv (J/mol*K)	
O_CP =  squeeze(OXYGEN(:,9,:)); %Cp (J/mol*K)	
O_SOUND_SPD =  squeeze(OXYGEN(:,10,:)); %Sound Spd. (m/s)	
O_JOULE_THOMPSON =  squeeze(OXYGEN(:,11,:));  %Joule-Thomson (K/MPa)	
O_VISCOSITY =  squeeze(OXYGEN(:,12,:));  %Viscosity (Pa*s)	
O_THERMAL_CONDUCTIVITY =  squeeze(OXYGEN(:,13,:)); %Therm. Cond. (W/m*K)	
O_PHASE =  squeeze(OXYGEN(:,14,:));%Phase

%%
%%Isolating Supercrtical Region
O_DENSITY(1:25,:) =  [];
O_VOLUME(1:25,:) =  [];	
O_INTERNAL_ENERGY(1:25,:) =  [];  
O_ENTHALPY(1:25,:) =  [];	
O_ENTROPY(1:25,:) =  [];
O_CV(1:25,:) =  [];
O_CP(1:25,:) =  [];
O_SOUND_SPD(1:25,:) =  []; 
O_JOULE_THOMPSON(1:25,:) =  [];
O_VISCOSITY(1:25,:) =  [];
O_THERMAL_CONDUCTIVITY(1:25,:) =  [];
O_PHASE(1:25,:) =  [];

O_DENSITY(379,:) =  [];
O_VOLUME(379,:) =  [];	
O_INTERNAL_ENERGY(379,:) =  [];  
O_ENTHALPY(379,:) =  [];	
O_ENTROPY(379,:) =  [];
O_CV(379,:) =  [];
O_CP(379,:) =  [];
O_SOUND_SPD(379,:) =  []; 
O_JOULE_THOMPSON(379,:) =  [];
O_VISCOSITY(379,:) =  [];
O_THERMAL_CONDUCTIVITY(379,:) =  [];
O_PHASE(379,:) =  [];

%%
%Extend H_density using cubic spline function

H_DENSITY(H_DENSITY==0) = NaN;
F = fillmissing(H_DENSITY,'pchip',2,'EndValues','extrap');
H_DENSITY = fillmissing(F,'pchip',1,'EndValues','extrap');

%importfile(H_GIBBS.mat);
%importfile(O_GIBBS.mat);

