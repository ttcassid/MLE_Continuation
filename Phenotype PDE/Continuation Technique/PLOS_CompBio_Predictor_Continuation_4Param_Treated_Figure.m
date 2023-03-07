%Script to simulate parameter sensitivity to data
% This script runs multiple steps of the predictor corrector algorithm but
% not continuation in the classic sense of solving the nonlinear system to
% find the optimal next point.

close all
clear all
rng(3) 

%% Initialize the optimizer
% tic
opts = optimset('MaxFunEvals',4000, 'MaxIter',5000,'Display','iter','TolFun', 5e-4,'TolX',1e-4);
X0 = [ 0.482725385475167, 0.349779478294077,0.419814106423369 , 0.702527043639005] ; %Parameters to be fit PA.Ra PA.Rb, PA.Da=PA.Db, PA.DaMax.
D =  [1616,8954,20864,52431,5683,6536] ;
TimeVec = [0,2,4,6].*24;
lb = [1e-2,1e-2 ,1e-2, 1e-2]; %Lower bounds
ub = [0.6, 0.48,0.48,1.05]; %Upper bounds
 
XInitialSol = X0; 
%% Set up perturbed data

% Successful run refitting d_a^treated
NPoints = 10 ; %  5 ;
PerturbNormLowerLimit = 0.01; % 0.05;
TotalPerturb = 0.75; % The total perturbation from the first data point.

PerturbStepSize =   (TotalPerturb-PerturbNormLowerLimit.*NPoints)*2./(NPoints*(NPoints+1)); 
PerturbNorm =  PerturbNormLowerLimit +PerturbStepSize.*[1:NPoints];  
DataPerturbationUnscaled = -1+2.*rand(NPoints,length(D)-1);
DataPerturbation  = zeros(NPoints,length(D)-1) ;
DataVec = zeros(NPoints+1,length(D));
DataVec(1,:) = log(D)./log(10);
for ii = 1:NPoints
    %Ensure that the treated data and untreated data move in same direction
    if DataPerturbationUnscaled(ii,2)*DataPerturbationUnscaled(ii,4) <0 ;
        DataPerturbationUnscaled(ii,2) = - DataPerturbationUnscaled(ii,2);
    end
    if DataPerturbationUnscaled(ii,2)*DataPerturbationUnscaled(ii,3) <0 ;
        DataPerturbationUnscaled(ii,3) = - DataPerturbationUnscaled(ii,3);
    end
    if DataPerturbationUnscaled(ii,3)*DataPerturbationUnscaled(ii,5) <0 ;
        DataPerturbationUnscaled(ii,5) = - DataPerturbationUnscaled(ii,5);
    end
    DataPerturbation(ii,:) = PerturbNorm(ii).*( DataPerturbationUnscaled(ii,:)./norm(DataPerturbationUnscaled(ii,:)) );
    DataVec(ii+1,:) = (ones(1,length(D)) + [0,DataPerturbation(ii,:)]).*DataVec(ii,:); 
end

%% Set up predictor
StepSize = 0.01;  
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

for jj = 1: NPoints
    %% Perturb data
    D1 =  (ones(1,length(D)) + [0,DataPerturbation(jj,:)]).*log(D)./log(10);   
    %% Calculate numerical gradient
    %Calculate mixed partial d/dD (d/dx G(x,D)) 
    for ii = 1 :NParam
        XLower(ii,:)    =  XUpdatedSol - StepSize.*ParameterStencil(ii,:); 
        XUpper(ii,:)    =  XUpdatedSol + StepSize.*ParameterStencil(ii,:); 
        FValLower(ii,:) =  PLOS_CompBio_Treated_ModelSolver_4Param(XLower(ii,:),D) ; 
        FValUpper(ii,:) = PLOS_CompBio_Treated_ModelSolver_4Param(XUpper(ii,:),D); 
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
                FValPositive(kk,:) = PLOS_CompBio_Treated_ObjFunction_4Param(XPositive(kk,:),D);
                FValNegative(kk,:) = PLOS_CompBio_Treated_ObjFunction_4Param(XNegative(kk,:),D);
            end
            ObjHessian(ii,mm) = ( FValPositive(1,:)+FValPositive(2,:) - ( FValNegative(1,:)+FValNegative(2,:) ) )./(4*StepSize.^2);
        end
    end

    % As we haven't calculated diagonal elements yet and enforce symmetry
    ObjHessian = ObjHessian + ObjHessian';
    
     %Calculate Diagonal elements of Hessian Matrix d^2/dx^2(G(x,D))
    for ii = 1:NParam
        ObjHessian(ii,ii) = (PLOS_CompBio_Treated_ObjFunction_4Param(XUpper(ii,:),D) - 2.*PLOS_CompBio_Treated_ObjFunction_4Param(XUpdatedSol,D) +PLOS_CompBio_Treated_ObjFunction_4Param(XLower(ii,:),D) )./(StepSize.^2);
    end 
    
    %Calculate the predictor matrix
    % dThetaDData = - inv(ObjHessian)*ObjMixedSecondDifferenceXD
%     cond(ObjHessian)
    dThetaDData = - ObjHessian\ObjMixedSecondDifferenceXD; 
    
    Corrector = dThetaDData*(D1'-log(D')./log(10))
    UpdatedGuess(jj,:) = XUpdatedSol + Corrector' ;
    % Update the data and objective function
    DUpdate = 10.^(D1);
    ObjFunction = @(X) PLOS_CompBio_Treated_ObjFunction_4Param(X,DUpdate);
    %Calculate  the Objective funciton at the updated guess and old
    %optimized point
    
    PredictedObjFunction(jj) = ObjFunction(UpdatedGuess(jj,:));
    GuessObjFunction(jj) = ObjFunction( XUpdatedSol );
    
    Test = PredictedObjFunction(jj)
    NaiveGuess = GuessObjFunction(jj)
    
    %Refit the updated data
    [XUpdate,fvalUpdate,~,output]=  fmincon(ObjFunction,XUpdatedSol,[],[],[],[],lb,ub,[],opts) ;
    XUpdateVec(jj,:) = XUpdate;
    FValUpdateVec(jj) = fvalUpdate;
    IterationsNaive(jj) = output.iterations;
    fEvalNaive(jj) = output.funcCount;
    
    [XCorrected,fvalCorrected,~,output] =  fmincon(ObjFunction,UpdatedGuess(jj,:),[],[],[],[],lb,ub,[],opts) ;
    
    XCorrectedVec(jj,:) = XCorrected;
    FValCorrectedVec(jj) = fvalCorrected;
    IterationsUpdate(jj) = output.iterations;
    fEvalUpdate(jj) = output.funcCount;
    % Reset the algorithm using the perturbed data point as the new starting point
    XUpdatedSol = XCorrected;
    
    D = DUpdate;
    PredictionfEvals(jj) = 2*NParam*(NParam+2);
    
end
TestVecAccuracy = PredictedObjFunction < GuessObjFunction;
TestVecIterationEffiency = IterationsUpdate < IterationsNaive;
TestVecFuncEval = fEvalUpdate < fEvalNaive;

TotalFEvalUpdate = sum(fEvalUpdate);
TotalFEvalNaive =  sum(fEvalNaive);

TotalIterUpdate = sum(IterationsUpdate);
TotalIterNaive =  sum(IterationsNaive);

PlottingPerturbNorm = zeros(1,length(PerturbNorm)) ;
PlottingIterUpdate = zeros(1,length(PerturbNorm)) ;
PlottingIterNaive = zeros(1,length(PerturbNorm)) ;
PlottingFEvalUpdate = zeros(1,length(PerturbNorm)) ;
PlottingFEvalNaive = zeros(1,length(PerturbNorm)) ;
PlottingFEvalPrediction = zeros(1,length(PerturbNorm)) ;

for ii = 1 :length(PerturbNorm)
    PlottingPerturbNorm(ii) = sum( PerturbNorm(1:ii));
    
    PlottingIterUpdate(ii) = sum( IterationsUpdate(1:ii));
    PlottingIterNaive(ii) = sum( IterationsNaive(1:ii));
    
    PlottingFEvalUpdate(ii) = sum( fEvalUpdate(1:ii));
    PlottingFEvalNaive(ii) = sum( fEvalNaive(1:ii));
    
    PlottingFEvalPrediction(ii) = sum( PredictionfEvals(1:ii));
    
end

LineWidth = 1.75;

figure(1)
g1 = plot(PlottingPerturbNorm(:),PredictedObjFunction(:) ,'--o','LineWidth',LineWidth,'Color',[67,147,195]/256 ,'MarkerSize',2,'MarkerFaceColor',[171,217,233]/256);
hold on
g2 = plot(PlottingPerturbNorm(:),GuessObjFunction(:) ,'--d','LineWidth',LineWidth,'Color',[239,138,98]/256,'MarkerSize',2,'MarkerFaceColor',[239,138,98]/256);
hold on
g4 = plot(PlottingPerturbNorm(:),FValCorrectedVec(:) ,'-x','LineWidth',LineWidth,'Color',[67,147,195]/256,'MarkerSize',2,'MarkerFaceColor',[213,62,79]/256);
hold on
ylim([0.1*min( PredictedObjFunction(:)) 1.2*max( GuessObjFunction(:) )])
xlim([0.8.*min(PlottingPerturbNorm(:))  1.02.*max(PlottingPerturbNorm(:))]);

figure(3)
g1 = plot([1:NPoints],PredictedObjFunction(:) ,'--o','LineWidth',LineWidth,'Color',[67,147,195]/256 ,'MarkerSize',2,'MarkerFaceColor',[171,217,233]/256);
hold on
g2 = plot([1:NPoints],GuessObjFunction(:) ,'--d','LineWidth',LineWidth,'Color',[239,138,98]/256,'MarkerSize',2,'MarkerFaceColor',[239,138,98]/256);
hold on
g4 = plot([1:NPoints],FValCorrectedVec(:) ,'-x','LineWidth',LineWidth,'Color',[67,147,195]/256,'MarkerSize',2,'MarkerFaceColor',[213,62,79]/256);
hold on
g3 = plot([1:NPoints],FValUpdateVec(:) ,'--x','LineWidth',LineWidth,'Color',[179,88,6]/256,'MarkerSize',2,'MarkerFaceColor',[213,62,79]/256);
hold on
ylim([0.1*min( PredictedObjFunction(:)) 1.2*max( GuessObjFunction(:) )])
xlim([0 length(NPoints)]);

figure(4)
g1 = plot([1:NPoints],PlottingIterUpdate(:) ,'-o','LineWidth',LineWidth,'Color',[67,147,195]/256 ,'MarkerSize',2,'MarkerFaceColor',[171,217,233]/256);
hold on
g2 = plot([1:NPoints],PlottingIterNaive(:) ,'-d','LineWidth',LineWidth,'Color',[239,138,98]/256,'MarkerSize',2,'MarkerFaceColor',[239,138,98]/256);

figure(5)
g1 = plot([1:NPoints],PlottingFEvalUpdate(:) ,'-o','LineWidth',LineWidth,'Color',[67,147,195]/256 ,'MarkerSize',2,'MarkerFaceColor',[171,217,233]/256);
hold on
g2 = plot([1:NPoints],PlottingFEvalNaive(:) ,'-d','LineWidth',LineWidth,'Color',[239,138,98]/256,'MarkerSize',2,'MarkerFaceColor',[239,138,98]/256);
hold on 
g3 = plot([1:NPoints],PlottingFEvalPrediction(:) ,'-s','LineWidth',LineWidth,'Color',[171,221,164]/256,'MarkerSize',2,'MarkerFaceColor',[171,221,164]/256);
hold on 



