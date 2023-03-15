function [ModelValues] =  Generic_ModelSimulation(X,D,T);

SampleTimeVec =  T ; 
%% Calculate model parameters


%% Time set up
ts = 0;
tf = SampleTimeVec(end); 
TotalTime = [ts tf];

%% Inital Conditions 
% Input the initial conditions

%% Solve the ODE systems

ModelValues =[]; % evaluate the solution of the ODE (following transformation, if necessary) at the sample times.
% Output Model Values.
 

%% ODE Solvers 


 


end