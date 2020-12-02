% script to demonstrate Lagrangian nature of the irradiance memory property
% in a one dimensional mixed layer.  Based on model of Kida & Ito (2017)

% started 15/10/2020 during response to reviewer comments.

% ---------------------%

%%%  Declare constants  %%%

% time step (s)
dt=60;
% diffusivity parameter (m^2/s)
%Kz=0.01; % FAST MIXING
Kz=0.001; %SLOW MIXING
% light attenuation coefficient (1/m)
k=0.04;
% mixed layer depth (m)
D=20;
% time scale to acclimate to irradiance (#timesteps)
gamma = 1440; % --- SLOW ACCLIMATION = 1day ---
%gamma = 60; % --- FAST ACCLIMATION = 1 hour ---
% surface irradiance (W/m^2)
I0=100;
% number of timesteps
niter=5000;
% number of phytoplankton cells
npart=100;

% ---------------------%

%%%  Declare variables %%%

% phytoplankton cell position
z(1:npart,1:niter)     = 0;
% instantaneous irradiance
I(1:npart,1:niter)     = 0;
% irradiance memory
I_mem(1:npart,1:niter) = 0;
% position memory (theoretical)
z_mem(1:npart,1:niter) = 0;

% ---------------------%

%%%  Run model %%%

% calculate independent trajectory for each phytoplankton cell
for i=1:npart

% initialise depth, using uniformly distributed random number
z(i,1) = D*rand;

    % loop over number of timesteps
    for t=2:niter
        % diffuse to new depth, using normally distributed random number
        z(i,t)=z(i,t-1)+normrnd(0,1)*(2*Kz*dt)^0.5; 
        % motion reflected at sea surface
        if (z(i,t)<0)
            z(i,t)=-z(i,t); 
        end
        % motion reflected at mixed layer base
        if (z(i,t)>D)
        z(i,t)=2*D-z(i,t); 
        end
        % calculate instantaneous irradiance
        I(i,t)=I0*exp(-k*z(i,t));
        if(t==2)
            % initialise irradiance memory
            I_mem(i,t)=I(i,t);
        else
            % update irradiance memory using acclimation timescale gamma
            I_mem(i,t)=I_mem(i,t-1)+(I(i,t)-I_mem(i,t-1))/gamma;
            % update position memory
            z_mem(i,t)=-log(I_mem(i,t)/I0)/k;
        end
    end

end
