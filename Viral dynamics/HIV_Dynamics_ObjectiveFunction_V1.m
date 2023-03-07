%% This is the objective function for the fitter for PA.Rb and PA.Da=PA.Db

function [Obj] = HIV_Dynamics_ObjectiveFunction_V1(X,D,T);

SampleTimeVec = T; 
%% Calculate model parameters
Y = exp(X);
PA.beta = Y(1);
PA.p = Y(2);
PA.delta = Y(3);
PA.c =  Y(4);
PA.N =  Y(5);
PA.lambda =  Y(6);

%% Time set up
ts = 0;
tf = SampleTimeVec(end); 
TotalTime = [ts tf];

%% Inital Conditions 
TIC = 180;
IIC = 20;
VIC = 50000;

%% Solve the ODE systems
ViralDynamicsIC = [TIC,IIC,VIC];

%% Solve the untreated system
 [sol1] = HIV_DynamicsODESolver(TotalTime,ViralDynamicsIC,PA); 
ModelValues = log( deval(sol1,SampleTimeVec,3))./log(10);
%% Calculate the Objective Function
Obj =  sqrt( sum((ModelValues-D).^2 ) )  ;

%% Solvers 
function [sol] = HIV_DynamicsODESolver(totaltime,IC,PA) %DDE model without therapy
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2);
sol = ode45(@HIV_DynamicsODE,totaltime,IC,opts);
function dydt = HIV_DynamicsODE(t,y); % y(1) = T(t), y(2) = I(t), y(3) = V(t);
dydt(1) = PA.lambda - PA.p*y(1) - PA.beta.*y(1).*y(3) ; %Differential equation for Abar(t)
dydt(2) = PA.beta.*y(1).*y(3) - PA.delta.*y(2) ; %Differential equation for Bbar(t)
dydt(3) = PA.delta.*PA.N.*y(2) - PA.c.*y(3) ;
dydt = dydt';
end

end

end