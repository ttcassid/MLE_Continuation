%Script to simulate parameter sensitivity to data for HIV dynamics model
close all
clear all
%rng default % For reproducibility
rng(1)

%% Initialize the optimizer
% tic
opts = optimset('MaxFunEvals',4000, 'MaxIter',5000,'Display','iter','TolFun', 1e-4,'TolX',1e-4);
%Parameters to be fit from viral dynamics model
X0 = log( [2.0e-05,0.15,0.55,5.5,900,80] ) ; %Known parameters taken from viral dynamics model
lb = log( [2.0e-06,0.05,0.1,2.0,500,40])  ; %Lower bounds
ub = log( [2.0e-04,0.5,1.0,8.0,1500,120])  ; %Upper bounds

% HIV Viral load data from proper simulation and perturbed with normally
% distributed noise.
D = [4.69897000433602,4.47550824516054,4.26977158768498,3.93827100436057,3.73611413009464,...
    3.63594058404862,3.58779862538249,3.49791944205020,3.41768427904502,3.40705458108154,...
    3.44292310115715,3.50618012621502,3.58064028389989,3.65294233543603,3.71322240212212,...
    3.75599631885288,3.78037311803958,3.78910543477391,3.78685161001207,3.77848369631436,...
    3.76800398294940,3.75813209009514,3.75034956546040,3.74516384237612,3.74242944421434,...
    3.74164123234798,3.74216426721666,3.74339168659599,3.74483561152213,3.74616304705651,...
    3.74719199297753,3.74786346452322,3.74820331057892,3.74828431910383,3.74819509161413,3.74801847186371,3.74781952314275] ;
TimeVec = [0,0.1,0.2:0.2:1,2:2:60] ;

% To copnsider only a subset of the simulated data.
DataToFitIndex = [4,7,11,14,17,25,30,36]; 

D = D(DataToFitIndex);
TimeVec = TimeVec(DataToFitIndex);

XInitialSol =  log( [2.0e-05,0.15,0.55 ,5.5,900,80] );
%% Set up perturbed data
NPoints = 4 ; 

SizePerturb = [0.1,0.2]; 
PerturbNormVec = [SizePerturb(1),-SizePerturb(1), SizePerturb(2),-SizePerturb(2)]; % linspace(PerturbNormLowerLimit,PerturbNormUpperLimit,NPoints);; % linspace(PerturbNormLowerLimit,PerturbNormUpperLimit,NPoints);
DataVec  = zeros(NPoints,length(D)) ;

%% Set up predictor
%Stepsize for the finite difference calculations
StepSize = 0.1; %0.1 Parameters are on the log scale so this is an approx ~10% fold change
NParam = length(X0);

%Structures necessary to calculate first order derivative
ParameterStencil = eye(NParam,NParam);
XLower = zeros(NParam,NParam);
XUpper = zeros(NParam,NParam);
FValLower = zeros(NParam,length(D));
FValUpper = zeros(NParam,length(D));
ObjMixedSecondDifferenceXD = zeros(NParam,length(D));

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

XCorrectedVec = zeros(NPoints,length(X0));
IterationsUpdate = zeros(NPoints,1);
fEvalUpdate = zeros(NPoints,1);
FValUpdate  = zeros(NPoints,1);
PredictionfEvals = zeros(NPoints,1);
ParameterAbsChangeVec =  zeros(NPoints,NParam);
BICDifferenceTableValue =  zeros(NPoints,1);

%% For plotting
Fig1 = figure(1);
Fig2 = figure(2);
%Colors
ColorVec = zeros(NPoints,3);
ColorVec(1,:) = [158,1,66]/256;
ColorVec(2,:) = [50,136,189 ]/256; % [213,62,79]/256;
ColorVec(3,:) = [244,109,67]/256;
ColorVec(4,:) = [94,79,162 ]/256; % [ 253,174,97]/256;
ColorVec(5,:) = [ 254,224,139]/256;
ColorVec(6,:) = [53,151,143]./256; % [ 230,245,152]/256;
ColorVec(7,:) = [ 171,221,164]/256;
ColorVec(8,:) = [ 102,194,165]/256;
ColorVec(9,:) = [50,136,189 ]/256;
ColorVec(10,:) = [94,79,162 ]/256;


for jj = 1: NPoints
    %% Perturb data 
     D1 =  D + PerturbNormVec(jj).*(abs( randn(1,length(D)) ) );   
    %% Calculate numerical gradient
    
    for ii = 1 :NParam
        XLower(ii,:)    =  XUpdatedSol - StepSize.*ParameterStencil(ii,:); 
        XUpper(ii,:)    =  XUpdatedSol + StepSize.*ParameterStencil(ii,:);
        FValLower(ii,:) =  HIV_Dynamics_ModelSolver_V1(XLower(ii,:),D,TimeVec) ;
        FValUpper(ii,:) = HIV_Dynamics_ModelSolver_V1(XUpper(ii,:),D,TimeVec);
        ObjMixedSecondDifferenceXD(ii,:) = -2.*(FValUpper(ii,:) - FValLower(ii,:))./(2.*StepSize);
    end
    
    % Create off diagonal elements of Hessian Matrix first
    %Calculate Diagonal elements of Hessian Matrix d^2/dx^2(G(x,D))
    %first, re-initalize to zero from previous runs
    ObjHessian = zeros(NParam,NParam);
    for ii = 1:NParam-1
        for mm = ii+1:NParam
            HessianIndex = FStencil([ii,mm]);
            XPositive(1,:) =  XUpdatedSol + StepSize.*PositiveStencil(HessianIndex,:);
            XPositive(2,:) =  XUpdatedSol + StepSize.*PositiveStencil(NEvalsHessian+HessianIndex,:);
            XNegative(1,:) =  XUpdatedSol + StepSize.*NegativeStencil(HessianIndex,:);
            XNegative(2,:) =  XUpdatedSol + StepSize.*NegativeStencil(NEvalsHessian+HessianIndex,:);
            for kk = 1:2
                FValPositive(kk,:) = HIV_Dynamics_ObjectiveFunction_V1(XPositive(kk,:),D,TimeVec);
                FValNegative(kk,:) = HIV_Dynamics_ObjectiveFunction_V1(XNegative(kk,:),D,TimeVec);
            end
            ObjHessian(ii,mm) = ( FValPositive(1,:)+FValPositive(2,:) - ( FValNegative(1,:)+FValNegative(2,:) ) )./(4*StepSize.^2);
        end
    end
    
    % As we haven't calculated diagonal elements yet and enforce symmetry
    ObjHessian = ObjHessian + ObjHessian';
    
    %Calculate Diagonal elements of Hessian Matrix d^2/dx^2(G(x,D))
    for ii = 1:NParam
        ObjHessian(ii,ii) = (HIV_Dynamics_ObjectiveFunction_V1(XUpper(ii,:),D,TimeVec) - 2.*HIV_Dynamics_ObjectiveFunction_V1(XUpdatedSol,D,TimeVec) +HIV_Dynamics_ObjectiveFunction_V1(XLower(ii,:),D,TimeVec) )./(StepSize.^2);
    end
    
    %Calculate the predictor matrix
    dThetaDData = - ObjHessian\ObjMixedSecondDifferenceXD;
    
    Corrector = dThetaDData*(D1'-D');
    UpdatedGuess(jj,:) = XUpdatedSol + Corrector' ;
    % Update the data and objective function
    DUpdate = D1;
    ObjFunction = @(X) HIV_Dynamics_ObjectiveFunction_V1(X,DUpdate,TimeVec);
    DataVec(jj,:) = D1;
    %Calculate  the Objective funciton at the updated guess and old
    %optimized point
    
    PredictedObjFunction(jj) = ObjFunction(UpdatedGuess(jj,:)); 
    
    Prediction = PredictedObjFunction(jj);
    XCorrectedVec(jj,:) = XUpdatedSol;
    
    [XCorrected,fvalCorrected,~,output] =  fmincon(ObjFunction,XUpdatedSol,[],[],[],[],lb,ub,[],opts) ;
    
    XCorrectedVec(jj,:) = XCorrected;
    IterationsUpdate(jj) = output.iterations;
    fEvalUpdate(jj) = output.funcCount;
    FValUpdate(jj) = fvalCorrected;
   PredictionfEvals(jj) = 2*NParam*(NParam+2);
%       
        %% Fold change in parameters
    ParameterAbsChangeVec(jj,:) = [( UpdatedGuess(jj,:) -X0 )] ; % [( UpdatedGuess(ii,:) -X0 ),norm(UpdatedGuess(ii,:) -X0) ];
    BICDifferenceTableValue(jj) = -2*Prediction+ 2*fvalCorrected ;%  -2*log(Prediction)+ 2*log(fvalCorrected) ;
end
disp('Logical test if the absolute difference in BIC is less than 2')
BICTest = abs(BICDifferenceTableValue) < 2

% PlottingPerturbNorm = zeros(1,length(PerturbNorm)) ;
PlottingIterUpdate = zeros(1,length(NPoints)) ; 
PlottingFEvalUpdate = zeros(1,length(NPoints)) ; 
PlottingFEvalPrediction = zeros(1,length(NPoints)) ;

for ii = 1 : NPoints
%     PlottingPerturbNorm(ii) = sum( PerturbNorm(1:ii));
    PlottingIterUpdate(ii) = sum( IterationsUpdate(1:ii));
    PlottingFEvalUpdate(ii) = sum( fEvalUpdate(1:ii)); 
    PlottingFEvalPrediction(ii) = sum( PredictionfEvals(1:ii));
end


%% Plotting scripts
PlottingIndex = 1:NPoints;
for nn = 1:  length(PlottingIndex);  
%     ii = 2*nn;
    ii = PlottingIndex(nn);
    %% Time set up
    ts = 0;
    tf = TimeVec(end);
    TotalTime = [ts tf];
    %Input model parameters from predicted guess
    Y = exp( UpdatedGuess(ii,:) );
    PA.beta = Y(1);
    PA.p = Y(2);
    PA.delta = Y(3);
    PA.c =   Y(4);
    PA.N =   Y(5);
    PA.lambda = Y(6);

    %% Inital Conditions
    TIC = 180;
    IIC = 20;
    VIC = 50000;
    %% Solve the ODE systems
    ViralDynamicsIC = [TIC,IIC,VIC];
    %% Solve the untreated system
    [solGuess] = HIV_DynamicsODESolver(TotalTime,ViralDynamicsIC,PA);
    %% Fig 1: Data vs. guessed solution
    figure(1)
    hold on
    g2 =  plot(solGuess.x(:),log(solGuess.y(3,:))./log(10),'LineWidth',1.1,'Color',ColorVec(ii,:),'LineStyle','-'); %grey
    hold on
    g12 = scatter(TimeVec,DataVec(ii,:),20,ColorVec(ii,:),'o','filled'); %,'20', [33,102,172]/255,'*');
    hold on
    ylim([2 5.5])
   %% Input model parameters from predicted guess
    Y = exp( XCorrectedVec(ii,:) );
    PA.beta = Y(1);
    PA.p = Y(2);
    PA.delta = Y(3);
    PA.c =   Y(4);
    PA.N =   Y(5);
    PA.lambda =   Y(6);
    %% Solve the untreated system
    [solFit] = HIV_DynamicsODESolver(TotalTime,ViralDynamicsIC,PA); 
    %%Fig 2: Data vs. fit solution
    figure(2)
    hold on
    g2 =  plot(solFit.x(:),log(solFit.y(3,:))./log(10),'LineWidth',1.1,'Color',ColorVec(ii,:),'LineStyle','-'); %grey
    hold on
    g12 = scatter(TimeVec,DataVec(ii,:),20,ColorVec(ii,:),'o','filled'); %,'20', [33,102,172]/255,'*');
    hold on
    ylim([2 5.5])
    
end

LineWidth = 1.75;
figure(3)
g1 = plot([1:NPoints],PlottingFEvalUpdate(:) ,'-o','LineWidth',LineWidth,'Color',[67,147,195]/256 ,'MarkerSize',3,'MarkerFaceColor',[67,147,195]/256 );
hold on
g3 = plot([1:NPoints],PlottingFEvalPrediction(:) ,'-s','LineWidth',LineWidth,'Color',[171,221,164]/256,'MarkerSize',3,'MarkerFaceColor',[171,221,164]/256);
hold on
xlim([0.95,4.05]);
xticks([1,2,3,4]);

 figure(4)
 h = heatmap(ParameterAbsChangeVec);

ObjFuncComparison =(PredictedObjFunction - FValUpdate')./FValUpdate'

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



