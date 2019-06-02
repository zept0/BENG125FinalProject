%%%%%%% BENG 125 Final Project 2019 %%%%%%%%
%%%%%%% Team Members %%%%%%%
% Jose Zaragoza, Allen Le, Jacqueline Olness, Kokila Perera, and Linh N Le %
% Date: 5/2019 %

%% Description of project:
% To model the trends of meme popularity using models defined in 
%Wang, L., & Wood, B. C. (2011). An epidemiological approach to model the viral propagation of memes. 
%Applied Mathematical Modelling, 35(11), 5442-5447. doi:10.1016/j.apm.2011.04.035

% Pseudo-code
% Preparing the model:
% 1. Define differential equations with unknown/unfitted parameters
%    Susceptible population: dS/dt = -alpha*S*I
%    Infected population: dI/dt = (alpha - beta)*S*I + (beta*N - gamma -
%    beta*I)*I
%    Parameters to tune: alpha = transmission rate prop cnst, beta = reinfection
%    rate prop cnst, gamma = loss of interest prop cnst
%    Constants: N = total population
% 2. Linearize the system using Jacobian (as dI/dt is non-linear in terms of I while system is in 2D)
% 3. Find eigenvalues and eigenvectors in terms of tuning parameters, I* &
%    S*
% 4. Load experimental data specific to a meme:
%        - population of infected individuals over time
%
% ***LOOP***
% 5. Define fixed points for the system: E = (S*, I*)
%       For a meme dying out: FP is (Sbar,0) where Sbar is a func of S0 &
%       I0
%       For a meme persisting: FP is (0,N - gamma/beta)
% 6. Solve eigenvalues & vectors for the equilibrium condition, E
% 7. Form equations for S and I as a function of time (using eigenvalues & eigenvectors)
% 8. Use initial conditions to solve for constants 
%    (this will be the final model equation where unknown parameters have to be tuned)
%
% Running modelling: 
% 9. Carry out least squares fitting on model equation with the
%    experimental data
% ***END LOOP***
%
% 10. Plot the infection curve with final tuned parameters
%    (eg. I vs t, or we can modify I as needed to represent some "search term index")
%%%%%%%%%%%%%%%%%%%%%%

clear all; 
clc;
close all; 
format long; 

%% Load experimental data


%% Run Model Fitting


%% Plot Results

