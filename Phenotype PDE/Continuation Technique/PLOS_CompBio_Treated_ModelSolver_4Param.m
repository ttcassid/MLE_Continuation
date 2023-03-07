function [ModelValues] =  PLOS_CompBio_Treated_ModelSolver_4Param(X,D);

UntreatedPopulationData = D(1:4); % WildTypePopulationData;
TreatedPopulationData = D([1,2,5,6]); 
 
WildTypeTime = [0,2,4,6].*24;
%% Calculate waning parameters SigmaA and SigmaB.
PA.PAmax = 0.95;  %Max probability of staying in phenotype A
PA.PAhomeo = 0.0;  %Long time probability of staying in phenotype A

PA.SigmaA = 1e-2; 
PA.SigmaB = 1e-2;

%% Parameters for treatment switching

%Taken from prior fitting and fixed
PA.Ra =  X(1);    %Reproductive rate A

PA.Rb = X(2);  %Reproductive rate B
PA.DaUntreated = X(3); %Clearance rate from population A
PA.DaTreated =  X(4);   % parameter to be fit- death rate in maximal drug concentration

PA.PBmax = 0.937275772950413; % X(5);  %Max probability of staying in phenotype B
PA.PBhomeo = 0.624850512178444 ; % X(6); %Long time probability of staying in phenotype B

PA.Db =  PA.DaUntreated;  % Clearance rate from population B


%% Carry capacity:
PA.CarryingCapacity = 1614384;
%% Relative Fitness Power
PA.n = 2; %2
PA.RelativeFitnessOnSwitch = 1; %Switch for relative fitness or not: 1 is on, 0 is off.
%% Time set up
ts = 0;
tf = 24.*7.5; 
TotalTime = [ts tf];

%% Inital Conditions
N = 10;
FLambda = @(lambda) 1- (4.*PA.Ra.*PA.Rb)./((PA.Ra+PA.DaUntreated+lambda).*(PA.Rb+PA.Db+lambda)) + ...
              ( 2.*PA.Rb./(PA.Rb+PA.Db+lambda)-1).*2.*PA.Ra.*( PA.PAhomeo./(PA.Ra+PA.DaUntreated+lambda)+(PA.PAmax-PA.PAhomeo)./(PA.Ra+PA.DaUntreated+PA.SigmaA+lambda) )+ ...
              ( 2.*PA.Ra./(PA.Ra+PA.DaUntreated+lambda)-1).*2.*PA.Rb.*( PA.PBhomeo./(PA.Rb+PA.Db+lambda)+(PA.PBmax-PA.PBhomeo)./(PA.Rb+PA.Db+PA.SigmaB+lambda) ) ; %Functional form for malthusian parameter x.
             options = optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
              LowerBound = max( min(PA.Rb-PA.Db,PA.Ra-PA.DaUntreated),-min(PA.Ra+PA.DaUntreated,PA.Rb+PA.Db)+0.001 );
              UpperBound =  max(PA.Rb-PA.Db,PA.Ra-PA.DaUntreated);
              StartingPoints = linspace(LowerBound,UpperBound,N);
         for nn = 1:N
            PossibleRoots(nn) = fzero(FLambda,StartingPoints(nn),options);
         end
Lambda = max(PossibleRoots);

SensitiveFrac = 0.9;
ANaught = SensitiveFrac.*UntreatedPopulationData(1).*(PA.Rb+PA.Db+Lambda);
ABarIC = SensitiveFrac.*UntreatedPopulationData(1);
BNaught =   (1-SensitiveFrac).*UntreatedPopulationData(1).*(PA.Rb+PA.Db+Lambda);
BBarIC = (1-SensitiveFrac).*UntreatedPopulationData(1);

NAAIC =  ANaught.*( PA.PAhomeo./(PA.Ra+PA.DaUntreated+Lambda)+(PA.PAmax-PA.PAhomeo)./(PA.Ra+PA.DaUntreated+PA.SigmaA+Lambda) );
NBBIC =  BNaught.*( PA.PBhomeo./(PA.Rb+PA.Db+Lambda)+(PA.PBmax-PA.PBhomeo)./(PA.Rb+PA.Db+PA.SigmaB+Lambda) );
TherapyIC = 0;
%% Solve the ODE systems
TreatCapacityIC = [ABarIC,BBarIC,NAAIC,NBBIC,TherapyIC];

%% Solve the untreated system
 [sol1] = AdaptiveSizeTimeDynamics(TotalTime,TreatCapacityIC,PA);
 Evalsol1Untreated = deval(sol1,WildTypeTime,1)+ deval(sol1,WildTypeTime,2);
 
 [sol1Treated] = AdaptiveSizeTreatedTimeDynamicsTreated(TotalTime,TreatCapacityIC,PA);
 EvalSol1Treated = deval(sol1Treated,WildTypeTime,1)+ deval(sol1Treated,WildTypeTime,2) ; %Calculate the number of cells at each time point.
 
ModelValues = [log( Evalsol1Untreated )./log(10), log( EvalSol1Treated([3,4]) )./log(10)];
 
%% Solvers
% Untreated Solver
function [sol] = AdaptiveSizeTimeDynamics(totaltime,IC,PA) %DDE model without therapy
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
sol = ode15s(@AdaptiveSizeTimeDynamics,totaltime,IC,opts);
function dydt = AdaptiveSizeTimeDynamics(t,y); % y(1) = Abar, y(2) = Bbar, y(3) = N_AA (A to A) ; y(4) = N_BB; y(5) = C(t)
dydt(1) = -( PA.Ra.*Logisitic([y(1),y(2)],t,PA)+ PA.DaUntreated  ).*y(1)  ...
          + 2.*PA.Ra.*Logisitic([y(1),y(2)],t,PA).*y(3) ...
          + 2.*RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA).*( y(2)-y(4) ) ; %Differential equation for Abar(t)
dydt(2) = - (RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA)+PA.Db  ).*y(2)  ...
          + 2.*PA.Ra.*Logisitic([y(1),y(2)],t,PA).*( y(1)- y(3) )...
          + 2.*RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA).*y(4) ; %Differential equation for Bbar(t)
dydt(3) = PA.PAmax.* ( 2.*PA.Ra.*Logisitic([y(1),y(2)],t,PA).*y(3)+ 2.*RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA).*( y(2)-y(4)) )...
            - ( PA.Ra.*Logisitic([y(1),y(2)],t,PA)+PA.DaUntreated  ) .*y(3)...
            -PA.SigmaA.*y(3)+PA.SigmaA.*PA.PAhomeo.*y(1);
dydt(4) = PA.PBmax.* ( 2.*PA.Ra.*Logisitic([y(1),y(2)],t,PA).*( y(1)- y(3) ) + 2.*RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA).*y(4) )...
           - ( RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA)+PA.Db ).*y(4)-PA.SigmaB.*y(4)+PA.SigmaB.*PA.PBhomeo.*y(2); 
dydt(5) =  0; %SizingDoseAdaptive(PA,[y(1),y(2)],t) - PA.TreatElim.*y(5);
dydt = dydt';
end

end
% Treated Solver
function [sol] = AdaptiveSizeTreatedTimeDynamicsTreated(totaltime,IC,PA) %DDE model without therapy
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
sol = ode15s(@AdaptiveSizeTreatedTimeDynamics,totaltime,IC,opts);
function dydt = AdaptiveSizeTreatedTimeDynamics(t,y); % y(1) = Abar, y(2) = Bbar, y(3) = N_AA (A to A) ; y(4) = N_BB; y(5) = C(t)
dydt(1) = -( PA.Ra.*Logisitic([y(1),y(2)],t,PA)+ DeathRate(t,PA) ).*y(1)  ...
          + 2.*PA.Ra.*Logisitic([y(1),y(2)],t,PA).*( 1- SwitchingProbability(t,PA) ).*y(3) ...
          + 2.*RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA).*( y(2)-y(4) ) ; %Differential equation for Abar(t)
dydt(2) = - (RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA)+ PA.Db   ).*y(2)  ...
          + 2.*PA.Ra.*Logisitic([y(1),y(2)],t,PA).*( y(1)- ( 1-SwitchingProbability(t,PA) ).*y(3) )...
          + 2.*RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA).*y(4) ; %Differential equation for Bbar(t)
dydt(3) = PA.PAmax.* ( 2.*PA.Ra.*Logisitic([y(1),y(2)],t,PA).*( 1- SwitchingProbability(t,PA) ).*y(3)+ 2.*RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA).*( y(2)-y(4)) )...
            - ( PA.Ra.*Logisitic([y(1),y(2)],t,PA)+DeathRate(t,PA) ) .*y(3)...
            -PA.SigmaA.*y(3)+PA.SigmaA.*PA.PAhomeo.*y(1);
dydt(4) = PA.PBmax.* ( 2.*PA.Ra.*Logisitic([y(1),y(2)],t,PA).*( y(1)- ( 1- SwitchingProbability(t,PA) ).*y(3) ) + 2.*RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA).*y(4) )...
           - ( RelativeFitnessB([y(1),y(2)],PA).*Logisitic([y(1),y(2)],t,PA)+PA.Db  ).*y(4)-PA.SigmaB.*y(4)+PA.SigmaB.*PA.PBhomeo.*y(2); 
dydt(5) =  0; % SizingDoseAdaptive(PA,[y(1),y(2)],t) - PA.TreatElim.*y(5);
dydt = dydt';
end

end


function y = DeathRate(t,PA)
    if t < 3.*24
        y = PA.DaUntreated;
    else
        y = PA.DaTreated;
    end        
end

function y = SwitchingProbability(t,PA)
    if t < 3.*24
        y = 0;
    else
        y = 0; %PA.SwitchingProbability;
    end        
end




function y = Logisitic(x,t,PA)
 y = 1 - (x(1)+x(2))./PA.CarryingCapacity;
end

function y = RelativeFitnessB(x,PA);
y = PA.Rb + PA.RelativeFitnessOnSwitch.*(PA.Ra-PA.Rb).*( x(2)^PA.n./(x(1)^PA.n+ x(2)^PA.n) );
end


end