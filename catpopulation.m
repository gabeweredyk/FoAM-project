% Proportion of scraps available to superpreadtors
S1 = 1;
% Proportion of scraps available to mesopredators
S0 = 1;
% Contribution of each succesful superpreadator-prey interaction to feral
% cat population
gamma = 1;
% Contribution of each succesful mesopredator-prey interaction to rat
% population
gamma_R = 1;
% Predation rate for Cat-Bird interactions
mu_B = 1;
% Predation rate for Cat-Rat interactions
mu_R = 1;
% Predation rate for Rat-Bird interactions
eta_B = 1;
% Per-capita rate of Feral, Rat & Bird populations
r_F = 1;
r_R = 1;
r_B = 1;
% Average lifespan of a stray cat
T_S = 1;
% Lifespan increase by TNR to indvidual cat
v_Plus = 1;
% Ratio between T_S and artificially shortened life span of stray cats
alpha = 1;
% Time it takes for an individual feral cat to undergo TNR
T_TNR = 1;
% x(t) = [ F(t); N(t); R(t); B(t) ]



dXdt = @(x) ...
[
  ( x(1)*( (gamma * (mu_B * x(4) + mu_R * x(3)))/(S1 + x(4) + x(3)) + r_F - (alpha/T_S) - (1/T_TNR)));
  ( x(1)/T_TNR - x(2)/(v_Plus * T_S) );
  ( x(3)*( (-mu_R*(x(1) + x(2)))/(S1 + x(4) + x(3)) + (gamma_R * eta_B * x(4))/(S0 + x(4)) + r_R));
  ( -x(4)*( (mu_B*(x(1) + x(2)))/(S1 + x(4) + x(3)) + (eta_B * x(3))/(S0 + x(4)) - r_B ) )
];


