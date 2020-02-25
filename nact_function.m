function [Potential] = nact_function(To, J, Rxn)
%Electrode Performance in Reversible Solid Oxide Fuel Cells
%Values derived from publication
%activation polarization potential estimated as a funciton of current
%density, temperature, and reaction completion

load('alpha.mat');
load('io.mat');
% To  Units [K]
% J   Units [A/cm^2];
% Rxn Units [0-1]

io_s   = io.fittedmodel(1/(To-273.15), Rxn);
alpha_s = alpha.fittedmodel(To -273.15,Rxn);

syms x
F = 96485.33289; %C mol?1
R = 8.314459;  %J?mol?1?K?1

%tofel equation
a1 = (1-alpha_s)*F/(R*To);
a2 = -alpha_s*F/(R*To);

eqn = io_s*(exp(a1*x)-exp(a2*x)) == J;

Potential = double(vpasolve(eqn, x));
end

