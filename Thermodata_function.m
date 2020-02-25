function [Cp_, H_, S_ ] = Thermodata_function(temperature,material)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

temperature = temperature + 273.15;
t = (temperature)/1000;

%http://webbook.nist.gov/cgi/cbook.cgi?Name=hydrogen&Units=SI&cTG=on&cTC=on&cTP=on#Thermo-Gas
%http://webbook.nist.gov/cgi/cbook.cgi?Name=water&Units=SI&cTG=on&cTC=on&cTP=on#Thermo-Gas
%http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas

%Oxygen

%100. - 700. Kelvin
%700. - 2000. Kelvin
%2000. - 6000. Kelvin

%A	31.32234	30.03235	20.91111
%B	-20.23531	8.772972	10.72071
%C	57.86644	-3.988133	-2.020498
%D	-36.50624	0.788313	0.146449
%E	-0.007374	-0.741599	9.245722
%F	-8.903471	-11.32468	5.337651
%G	246.7945	236.1663	237.6185
%H	0.0	        0.0	        0.0

Oxygen = [31.32234 -20.23531 57.86644  -36.50624 -0.007374 -8.903471 246.7945 0.0;
          30.03235  8.772972 -3.988133 0.788313  -0.741599 -11.32468 236.1663 0.0;
          20.91111 10.72071  -2.020498 0.146449   9.245722  5.337651 237.6185 0.0];

%Hydrogen

%298. - 1000.	Kelvin
%1000. - 2500.	Kelvin
%2500. - 6000.  Kelvin

%A	33.066178	18.563083	43.413560
%B	-11.363417	12.257357	-4.293079
%C	11.432816	-2.859786	1.272428
%D	-2.772874	0.268238	-0.096876
%E	-0.158558	1.977990	-20.533862
%F	-9.980797	-1.147438	-38.515158
%G	172.707974	156.288133	162.081354
%H	0.0	0.0	0.0

Hydrogen = [33.066178 -11.363417 11.432816 -2.772874  -0.158558  -9.980797 172.707974 0.0;
            18.563083  12.257357 -2.859786  0.268238   1.977990  -1.147438 156.288133 0.0;
            43.413560  -4.293079  1.272428 -0.096876 -20.533862 -38.515158 162.081354 0.0];

%Water

%500. - 1700.  Kelvin
%1700. - 6000. Kelvin
%A	30.09200	41.96426
%B	6.832514	8.622053
%C	6.793435	-1.499780
%D	-2.534480	0.098119
%E	0.082139	-11.15764
%F	-250.8810	-272.1797
%G	223.3967	219.7809
%H	-241.8264	-241.8264

Water = [ 0.0     0.0       0.0       0.0       0.0         0.0      0.0       0.0;
         30.09200 6.832514  6.793435 -2.534480  0.082139 -250.8810 223.3967 -241.8264;
         41.96426 8.622053 -1.499780  0.098119 -11.15764 -272.1797 219.7809 -241.8264];
     
data = [ 0.0     0.0       0.0       0.0       0.0         0.0      0.0       0.0];

i = 0;

switch material
   case 'oxygen'
      data = Oxygen;
      if(temperature < 700)
          i = 1;
      elseif(temperature <2000)
          i = 2;
      else 
          i = 3;
      end
   case 'hydrogen'
      data = Hydrogen;
      if(temperature < 1000)
          i = 1;
      elseif(temperature <2500)
          i = 2;
      else 
          i = 3;
      end
   case 'water'
      data = Water;
end


Cp_ = data(i, 1) + data(i, 2)*t + data(i, 3)*power(t,2) + data(i, 4)*power(t,3)  + data(i, 5)*power(t,-2) ; 

H_ = H_298 +  data(i, 1)*t + data(i, 2)*power(t,2)/2 + data(i, 3)*power(t,3)/3 + data(i, 4)*power(t,4)/4 - data(i, 5)/t + data(i, 6) - data(i, 8);

S_ = data(i, 1)*ln(t) + data(i, 2)*t + data(i, 3)*power(t,2) + data(i, 4)*power(t,3)/3 - data(i, 5)*power(t,-2)/2 + data(i, 7);



end

