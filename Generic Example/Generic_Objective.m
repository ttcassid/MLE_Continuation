%% This is the objective function for the fitter for PA.Rb and PA.Da=PA.Db

function [Obj] = Generic_Objective(X,D,T); % Input parameters X, data D, and time of data collection T.

SampleTimeVec = T; 

%% Calculate model parameters

%% Time set up
ts = 0;
tf = SampleTimeVec(end); 
TotalTime = [ts tf];

%% Inital Conditions 
% Set initial conditions for the ODE system
%% Solve the ODE system
%Solve the ODE system

ModelValues =; % evaluate the solution of the ODE (following transformation, if necessary) at the sample times.
%% Calculate the Objective Function
Obj =  sqrt( sum((ModelValues-D).^2 ) )  ;

%% ODE Solvers 


end