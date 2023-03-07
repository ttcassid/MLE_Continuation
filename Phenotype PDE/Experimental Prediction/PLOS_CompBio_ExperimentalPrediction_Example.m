%Script to simulate parameter sensitivity to data
% This script runs multiple steps of the predictor corrector algorithm but
% not continuation in the classic sense of solving the nonlinear system to
% find the optimal next point.

close all
clear all

%rng default % For reproducibility
rng(4)

%% Initialize the optimizer
% tic
opts = optimset('MaxFunEvals',4000, 'MaxIter',5000,'Display','iter','TolFun', 1e-5,'TolX',1e-4);
X0 = [ 0.482725385475167, 0.349779478294077,0.419814106423369 , 0.702527043639005] ; %Parameters to be fit PA.Ra PA.Rb, PA.Da=PA.Db, PA.DaMax.
ExperimentalData = [1616,8954,20864,52431,5683,6536] ; 
% Simulated data taken from model simulation with fixed parameters.
% Vector is repeated to account for perturbations by plus/minus the perturbation factor 
SimulatedData =  [8129.72,6344.82,5484.14, 5091.84, 4927.49, 5579.04, 6890.37, 8129.72,6344.82,5484.14, 5091.84, 4927.49, 5579.04, 6890.37 ]; 
D =  ExperimentalData ;
% [ExperimentalData(1:4),SimulatedData(1:5),ExperimentalData(5),SimulatedData(6), ExperimentalData(6),SimulatedData(end)];
 

ExperimentalTimeVec = [0,2,4,6].*24;
% Simulated measurement time taken from model simulation with fixed parameters.
% Vector is repeated to account for perturbations by plus/minus the perturbation factor 
AdditionalTimeVec = [3.1,3.2,3.3,3.4,3.5,5,7,3.1,3.2,3.3,3.4,3.5,5,7].*24;

TimeVec = ExperimentalTimeVec;  

 
lb = [1e-2, 1e-2,1e-2 ,1e-2]; %Lower bounds
ub = [0.6, 0.48,0.48,1.05]; %Upper bounds
XInitialSol = X0 ;  
%% Set up perturbed data 
NPoints = length(SimulatedData) ; 
PerturbStepSize =  0.3;  
DataPerturbation  = [SimulatedData(1:length(SimulatedData)/2)+PerturbStepSize.*8129.72,SimulatedData(1+length(SimulatedData)/2:end)-PerturbStepSize.*8129.72]; ; %[(1+PerturbStepSize).*SimulatedData(1:length(SimulatedData)/2),(1-PerturbStepSize).*SimulatedData(1+length(SimulatedData)/2:end)]; % perturb each simulated data point in positive/negative direction by the same step size
DataVec = zeros(NPoints+1,length(D)+1);
DataVec(1,:) = log([D,SimulatedData(1)])./log(10);

%% Set up predictor
StepSize = 0.01;  
NParam = length(X0);
%Structures necessary to calculate first order derivative
ParameterStencil = eye(NParam,NParam);
XLower = zeros(NParam,NParam);
XUpper = zeros(NParam,NParam);
% because we are appending one extra measurement onto each FVal evaluation,
% the length is |D|+1
FValLower = zeros(NParam,length(D)+1);
FValUpper = zeros(NParam,length(D)+1);
ObjMixedSecondDifferenceXD = zeros(NParam,length(D)+1);

% Structures necessary to calculate Hessian matrix
BasisVec = eye(NParam);
k = 1; % Counter for each row.
NEvalsHessian = (NParam)*(NParam-1)/2; 
FValPosTemplate = zeros(NEvalsHessian,NParam);
FValNegTemplate = zeros(NEvalsHessian,NParam);
for ii = 1 :NParam-1
    for jj = ii+1:NParam
        FValPosTemplate(k,:) = BasisVec(ii,:)+BasisVec(jj,:);
        FValNegTemplate(k,:) = -BasisVec(ii,:)+BasisVec(jj,:);
        k = k+1;
    end
end
PositiveStencil = [FValPosTemplate;-FValPosTemplate];
NegativeStencil = [FValNegTemplate;-FValNegTemplate];

% Function to identify the Stencil row for the i,j partial derivative;
% x(1)<x(2) necessarily here, only calculating the lower portion of the Hessian due to symmetry
FStencil = @(x) (x(1)-1)*(2*NParam - x(1)-2)./2 + x(2)-1;
% FStencilV1 = @(x) (x(1)-1)*(2*NParam - x(1))./2 + x(2);

%Matrices to store parameter updates and function evaluations
XPositive = zeros(NParam,NParam);
XNegative = zeros(NParam,NParam);
FValPositive = zeros(NParam,1);
FValNegative = zeros(NParam,1);
ObjHessian = zeros(NParam,NParam);

%% Set up output
UpdatedGuess = zeros(NPoints,length(X0));
XUpdatedSol =   XInitialSol;
PredictedObjFunction = zeros(NPoints,1);
GuessObjFunction = zeros(NPoints,1);
XUpdateVec = zeros(NPoints,length(X0));
FValUpdateVec = zeros(NPoints,1);
XCorrectedVec = zeros(NPoints,length(X0));
FValCorrectedVec = zeros(NPoints,1);

IterationsNaive = zeros(NPoints,1);
IterationsUpdate = zeros(NPoints,1);
fEvalNaive = zeros(NPoints,1);
fEvalUpdate = zeros(NPoints,1);

PredictionfEvals = zeros(NPoints,1);

ParameterFoldChangeVec = zeros(NPoints,length(X0)); % zeros(NPoints,length(X0)+1);
ParameterAbsChangeVec = zeros(NPoints,length(X0)); % zeros(NPoints,length(X0)+1);

% For plotting
Fig1 = figure(1);
Fig2 = figure(2);
%Colors
ColorVec = zeros(NPoints,3);
ColorVec(1,:) = [158,1,66]/256; 
ColorVec(2,:) = [213,62,79]/256;
ColorVec(3,:) = [244,109,67]/256;
ColorVec(4,:) = [158,1,66]/256; ; %[ 253,174,97]/256;
ColorVec(5,:) = [158,1,66]/256; ; % [ 254,224,139]/256;
ColorVec(6,:) = [53,151,143]./256; % [ 230,245,152]/256;
ColorVec(7,:) = [ 171,221,164]/256;
ColorVec(8,:) = [ 102,194,165]/256;
ColorVec(9,:) = [50,136,189 ]/256;
ColorVec(10,:) = [171,221,164]/256 ; 
ColorVec(11,:) = [171,221,164]/256 ;
ColorVec(12,:) = [171,221,164]/256;
ColorVec(13,:) = [171,221,164]/256;
ColorVec(14,:) = [171,221,164]/256;

for jj = 1:  NPoints
    %% Perturb data 
    D =  [ExperimentalData,SimulatedData(jj)];
    TimeVec = [ExperimentalTimeVec, AdditionalTimeVec(jj)];
    D1 = log([ExperimentalData,DataPerturbation(jj)])./log(10);


    %% Calculate numerical gradient
    %Calculate mixed partial d/dD (d/dx G(x,D))
   for ii = 1 :NParam
        XLower(ii,:)    =  XUpdatedSol - StepSize.*ParameterStencil(ii,:); 
        XUpper(ii,:)    =  XUpdatedSol + StepSize.*ParameterStencil(ii,:); 
        FValLower(ii,:) =  PLOS_CompBio_ExperimentalPredictionModelSolver(XLower(ii,:),D,TimeVec) ; 
        FValUpper(ii,:) = PLOS_CompBio_ExperimentalPredictionModelSolver(XUpper(ii,:),D,TimeVec); 
        ObjMixedSecondDifferenceXD(ii,:) = -2.*(FValUpper(ii,:) - FValLower(ii,:))./(2.*StepSize);
    end

   % Create off diagonal elements of Hessian Matrix first
    %Calculate Diagonal elements of Hessian Matrix d^2/dx^2(G(x,D)) 
    %first, re-initalize to zero from previous runs
    ObjHessian = zeros(NParam,NParam);
    for ii = 1:NParam-1
        for mm = ii+1:NParam
            HessianIndex = FStencil([ii,mm]);
            XPositive(1,:) =  XUpdatedSol + StepSize.*PositiveStencil(HessianIndex,:); %.*XInitialSol;
            XPositive(2,:) =  XUpdatedSol + StepSize.*PositiveStencil(NEvalsHessian+HessianIndex,:); %.*XInitialSol;
            XNegative(1,:) =  XUpdatedSol + StepSize.*NegativeStencil(HessianIndex,:); %.*XInitialSol;
            XNegative(2,:) =  XUpdatedSol + StepSize.*NegativeStencil(NEvalsHessian+HessianIndex,:); %.*XInitialSol;
            for kk = 1:2
                FValPositive(kk,:) = PLOS_CompBio_ExperimentalPrediction_ObjectiveFunction(XPositive(kk,:),D,TimeVec);
                FValNegative(kk,:) = PLOS_CompBio_ExperimentalPrediction_ObjectiveFunction(XNegative(kk,:),D,TimeVec);
            end
            ObjHessian(ii,mm) = ( FValPositive(1,:)+FValPositive(2,:) - ( FValNegative(1,:)+FValNegative(2,:) ) )./(4*StepSize.^2);
        end
    end

    % As we haven't calculated diagonal elements yet and enforce symmetry
    ObjHessian = ObjHessian + ObjHessian';
    
     %Calculate Diagonal elements of Hessian Matrix d^2/dx^2(G(x,D))
    for ii = 1:NParam
        ObjHessian(ii,ii) = (PLOS_CompBio_ExperimentalPrediction_ObjectiveFunction(XUpper(ii,:),D,TimeVec) ...
             - 2.*PLOS_CompBio_ExperimentalPrediction_ObjectiveFunction(XUpdatedSol,D,TimeVec)...
             +PLOS_CompBio_ExperimentalPrediction_ObjectiveFunction(XLower(ii,:),D,TimeVec) )./(StepSize.^2);
    end 
    
    %Calculate the predictor matrix
    dThetaDData = - ObjHessian\ObjMixedSecondDifferenceXD; 
    
    Corrector = dThetaDData*(D1'-log(D')./log(10))
    UpdatedGuess(jj,:) = XUpdatedSol + Corrector' ;
    % Update the data and objective function
    DUpdate = 10.^(D1);
    DataVec(jj,:) = DUpdate; 
    
    ObjFunction = @(X) PLOS_CompBio_ExperimentalPrediction_ObjectiveFunction(X,DUpdate,TimeVec);
    %Calculate  the Objective funciton at the updated guess and old
    %optimized point
    
    PredictedObjFunction(jj) = ObjFunction(UpdatedGuess(jj,:));
    GuessObjFunction(jj) = ObjFunction( XUpdatedSol );

end

%% Plotting

% plot the model trajectories
PlottingIndex = [4,NPoints/2+4]; 
for nn = 1: length(PlottingIndex);  
    
    ii =  PlottingIndex(nn);
     %% Time set up
    ts = 0;
    tf = 24.*7.5;
    TotalTime = [ts tf]; 
    WildTypeTime = [ExperimentalTimeVec, AdditionalTimeVec(ii)]  ;
    %% Data
    DUpdate = log([ExperimentalData,DataPerturbation(ii)])./log(10); 
    UntreatedPopulationData = 10.^( DUpdate(1:4) ) ; % WildTypePopulationData;
    TreatedPopulationData = 10.^( DUpdate([1,2,5:end]));
    
    %% Fixed parameters
    PA.PAmax = 0.95;  %Max probability of staying in phenotype A
    PA.PAhomeo = 0.0;  %Long time probability of staying in phenotype A
    PA.PBmax = 0.937275772950413; %%Max probability of staying in phenotype B
    PA.PBhomeo = 0.624850512178444 ; %  %Long time probability of staying in phenotype B
    
    PA.SigmaA = 1e-2;
    PA.SigmaB = 1e-2; 
    PA.CarryingCapacity = 1614384; 
    
    PA.n = 2; %2
    PA.RelativeFitnessOnSwitch = 1; %Switch for relative fitness or not: 1 is on, 0 is off.
    %% Parameters from sensitivity prediction
    X = UpdatedGuess(ii,:);
    PA.Ra = X(1) ; %  %Reproductive rate A
    PA.Rb = X(2);  %Reproductive rate B
    PA.DaUntreated = X(3); %Clearance rate from population A
    PA.DaTreated = X(4); % % parameter to be fit- death rate in maximal drug concentration   
    PA.Db =  PA.DaUntreated;  % Clearance rate from population B
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
    
    %Initial conditions
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
    % Solve the treated system
    [sol1Treated] = AdaptiveSizeTreatedTimeDynamicsTreated(TotalTime,TreatCapacityIC,PA);
    
    figure(10)
    hold on
     g1 =  plot(sol1.x(:)./24,   sol1.y(1,:)+sol1.y(2,:) ,'LineWidth',1.1,'Color',ColorVec(ii,:),'LineStyle','-'); %grey
    hold on
    g2 =  plot(sol1Treated.x(:)./24, sol1Treated.y(1,:)+sol1Treated.y(2,:) ,'LineWidth',1.1,'Color',ColorVec(ii,:),'LineStyle','--'); %grey
    hold on
    g12 = scatter(WildTypeTime(1:4)./24,  UntreatedPopulationData(:) ,40,'b','o','filled'); %,'20', [33,102,172]/255,'*');  
    hold on
     g22 = scatter(WildTypeTime([3,4])./24, TreatedPopulationData([3,4])  ,40,'b','o','filled') ; %,'20', [33,102,172]/255,'*');
    hold on
    TreatedPopulationData = 10.^( DUpdate([1,2,5:end]));
end

   % Plot the perturbed data points
   for ii = 4: 4;  
       WildTypeTime = [ExperimentalTimeVec, AdditionalTimeVec(ii)]  ;
       DUpdate = log([ExperimentalData,DataPerturbation(ii)])./log(10);
       TreatedPopulationData = 10.^( DUpdate([1,2,5:end]));
       hold on
       
       figure(10);
       TreatedPopulationData = 10.^( DUpdate([1,2,5:end]));
       g23 = scatter(WildTypeTime(5:end)./24, TreatedPopulationData(5:end)  ,40,ColorVec(ii,:),'d','filled') ; %,'20', [33,102,172]/255,'*');
       hold on
       % Plot both perturbations for figure
       DUpdateNegativePerturbation = log([ExperimentalData,DataPerturbation(NPoints/2 + ii)])./log(10);
       TreatedPopulationData = 10.^( DUpdateNegativePerturbation([1,2,5:end]));
       g24 = scatter(WildTypeTime(5:end)./24, TreatedPopulationData(5:end)  ,40,ColorVec(floor(NPoints/2)+(ii-1),:),'s','filled') ; %,'20', [33,102,172]/255,'*');
       hold on
   end
   

    %% Fold change in parameters
    for ii = 1:NPoints
        ParameterFoldChangeVec(ii,:) = [( UpdatedGuess(ii,:) -X0 )./X0 ];       
    end
    
    XSpace = [1:NPoints/2];
    TornadoValues = 100.*[ParameterFoldChangeVec(1:NPoints/2,4),ParameterFoldChangeVec(1+NPoints/2:end,4)];
    figure(4)
    TornadoColorVec1 =  [158,1,66]/256; 
    TornadoColorVec2 =  [ 171,221,164]/256;
   b =  barh(XSpace,TornadoValues,'stacked');
   b(1).FaceColor = TornadoColorVec1;
   b(2).FaceColor = TornadoColorVec2;
   xlim([-15 17])
   yticklabels([]); 
    
    figure(2)
    h = heatmap(ParameterFoldChangeVec);
    

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
dydt(5) =  0;
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


function y = Logisitic(x,t,PA)
 y = 1 - (x(1)+x(2))./PA.CarryingCapacity;
end

function y = RelativeFitnessB(x,PA);
y = PA.Rb + PA.RelativeFitnessOnSwitch.*(PA.Ra-PA.Rb).*( x(2)^PA.n./(x(1)^PA.n+ x(2)^PA.n) );
end


function y = SwitchingProbability(t,PA)
    if t < 3.*24
        y = 0;
    else
        y = 0; 
    end        
end


