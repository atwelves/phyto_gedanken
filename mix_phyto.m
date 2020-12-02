% Script to investigate apparent paradox of increased surface productivity
% under decreased light penetration

% Set number of iterations
max_iter = 200;

% Set maximum depth
depth = 20; %m

% Constants for phytoplankton growth, 
mu_max = 1.47; % maximum (nutrient and light saturated) per capita growth rate (1/day)
m = 0.19; % per capita mortality rate (1/day)
K_Fe = 1.6e-7; % half saturation constant for iron limitation (molC/m^3)
Irrk = 70; % photosynthesis efficiency constant - depends on iron in BLING
%theta = 0.01; % chlorophyll to carbon ratio

% Constants for nutrient (iron) supply
FetoN_min = 1.3e-5; % iron uptake ratio (mol Fe/mol C)
FetoN_max = 1.6875e-4
d_FetoN   = FetoN_max-FetoN_min;
k_FetoN   = 8e-10;
CtoN      = 81;
r = 0; % Per capita iron loss due to sinking etc. (1/day)
f = 0; % Iron flux into domain (mol Fe/m^3/day)

% Optical parameters
kappa = 0.05;
I0    = 100;

% Set initial conditions
P        (1:depth,1:max_iter) = 1e-10;
Fe       (1:depth,1:max_iter) = 1e-7;
FetoN    (1:depth,1:max_iter) = 0;
I        (1:depth,1:max_iter) = 0;
light_lim(1:depth,1:max_iter) = 0;
iron_lim (1:depth,1:max_iter) = 0;
mu       (1:depth,1:max_iter) = 0;
upt      (1:depth,1:max_iter) = 0;

for t = 1:max_iter
    for k = 1:depth
	if (Fe(k,t)<0)
		Fe(k,t) = 0;
    end
        % calculate local irradiance level
        I(k,t)         = I0*exp(-kappa*k);
        % calculate light and iron limitation terms
        light_lim(k,t) = 1-exp(-I(k,t)/Irrk);
        iron_lim (k,t) = Fe(k,t)/(Fe(k,t) + K_Fe);
        % calculate per capita growth rate
        mu(k,t)        = mu_max*light_lim(k,t)*iron_lim(k,t); 
        NPP(k,t)       = mu(k,t)*P(k,t);
        % calculate iron uptake
        upt_lim(k,t)   = Fe(k,t)/(Fe(k,t) + k_FetoN);
        FetoN(k,t)     = FetoN_min + d_FetoN*upt_lim(k,t);
        upt(k,t)       = FetoN(k,t)*mu(k,t)*P(k,t)/CtoN;
       
        % update iron and biomass vectors
        Fe(k,t+1)         = f + (1-r)*Fe(k,t) - upt(k,t); 
        P (k,t+1)         = P(k,t)*(1 - m + mu(k,t));
    end
    % OPTION: Homogenise iron by assuming that uptake is applied equally 
    %           to entire water column
    Fe(:,t+1)=mean(Fe(:,t+1));
    P (:,t+1)=mean(P (:,t+1));
end


