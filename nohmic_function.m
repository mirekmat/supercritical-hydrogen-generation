function [Potential] = nohmic_function(J,thickness_YSZ,T)
%NOHMIC_FUNCTION.m
%thickness (um)
Potential = (2.99*10^-7)*J*thickness_YSZ*exp(10300/T);
%volts
end

