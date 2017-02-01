function [resultGUI,info] = matRad_fluenceOptimization(dij,cst,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad inverse planning wrapper function
% 
% call
%   [resultGUI,info] = matRad_fluenceOptimization(dij,cst,pln)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   info:       struct containing information about optimization
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% issue warning if biological optimization impossible
if sum(strcmp(pln.bioOptimization,{'LEMIV_effect','LEMIV_RBExD','LSM_effect','LSM_RBExD'}))>0 && (~isfield(dij,'mAlphaDose') || ~isfield(dij,'mSqrtBetaDose')) && strcmp(pln.radiationMode,'carbon')
    warndlg('Alpha and beta matrices for effect based and RBE optimization not available - physical optimization is carried out instead.');
    pln.bioOptimization = 'none';
end

if ~isdeployed % only if _not_ running as standalone
    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'optimization'))
    
    % get handle to Matlab command window
    mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
    cw          = mde.getClient('Command Window');
    xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
    h_cw        = handle(xCmdWndView,'CallbackProperties');

    % set Key Pressed Callback of Matlab command window
    set(h_cw, 'KeyPressedCallback', @matRad_CWKeyPressedCallback);

end

% initialize global variables for optimizer
global matRad_global_x;
global matRad_global_d;
global matRad_STRG_C_Pressed;
global matRad_objective_function_value;
global matRad_iteration;
global kDVH;
global kDCH;
global JACOBIAN;
global GRADIENT;
global CONSTRAINT
global fScaling;
global cScaling;

matRad_global_x                 = NaN * ones(dij.totalNumOfBixels,1);
matRad_global_d                 = NaN * ones(dij.numOfVoxels,1);
matRad_STRG_C_Pressed           = false;
matRad_objective_function_value = [];
matRad_iteration                = 0;
kDVH                            = [];
kDCH                            = [];
JACOBIAN                        = [];
GRADIENT                        = [];
CONSTRAINT                      = [];
fScaling                        = 1;
cScaling                        = 1;

% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

% adjust objectives and constraints internally for fractionation 
for i = 1:size(cst,1)
    for j = 1:size(cst{i,6},1)
       cst{i,6}(j).dose = cst{i,6}(j).dose/pln.numOfFractions;
    end
end

% find target indices and described dose(s) for weight vector
% initialization
V          = [];
doseTarget = [];
ixTarget   = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;cst{i,4}{1}];
        doseTarget = [doseTarget cst{i,6}.dose];
        ixTarget   = [ixTarget i*ones(1,length([cst{i,6}.dose]))];
    end
end
[doseTarget,i] = max(doseTarget);
ixTarget       = ixTarget(i);
wOnes          = ones(dij.totalNumOfBixels,1);

% set the IPOPT options.
matRad_ipoptOptions;

% modified settings for photon dao
if pln.runDAO && strcmp(pln.radiationMode,'photons')
%    options.ipopt.max_iter = 50;
%    options.ipopt.acceptable_obj_change_tol     = 7e-3; % (Acc6), Solved To Acceptable Level if (Acc1),...,(Acc6) fullfiled

end

% set bounds on optimization variables
options.lb              = zeros(1,dij.totalNumOfBixels);        % Lower bound on the variables.
options.ub              = inf * ones(1,dij.totalNumOfBixels);   % Upper bound on the variables.
funcs.iterfunc          = @(iter,objective,paramter) matRad_IpoptIterFunc(iter,objective,paramter,options.ipopt.max_iter);
    
% calculate initial beam intensities wInit
if  strcmp(pln.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
    % check if a constant RBE is defined - if not use 1.1
    if ~isfield(dij,'RBE')
        dij.RBE = 1.1;
    end
    bixelWeight =  (doseTarget)/(dij.RBE * mean(dij.physicalDose{1}(V,:)*wOnes)); 
    wInit       = wOnes * bixelWeight;
        
elseif ((strcmp(pln.bioOptimization,'LEMIV_effect') || strcmp(pln.bioOptimization,'LEMIV_RBExD')) && strcmp(pln.radiationMode,'carbon')) || ...
       ((isequal(pln.bioOptimization,'LSM_effect')  || isequal(pln.bioOptimization,'LSM_RBExD'))  && strcmp(pln.radiationMode,'protons'))

    % check if you are running a supported rad
    dij.ax   = zeros(dij.numOfVoxels,1);
    dij.bx   = zeros(dij.numOfVoxels,1);
    
    for i = 1:size(cst,1)
        
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
             dij.ax(cst{i,4}{1}) = cst{i,5}.alphaX;
             dij.bx(cst{i,4}{1}) = cst{i,5}.betaX;
        end
        
        for j = 1:size(cst{i,6},2)
            % check if prescribed doses are in a valid domain
            if cst{i,6}(j).dose > 5 && isequal(cst{i,3},'TARGET')
                error('Reference dose > 5Gy[RBE] for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.');
            end
            
        end
    end
     
    if isequal(pln.bioOptimization,'LSM_effect') || isequal(pln.bioOptimization,'LEMIV_effect')
        
           effectTarget = cst{ixTarget,5}.alphaX * doseTarget + cst{ixTarget,5}.betaX * doseTarget^2;
           p            = (sum(dij.mAlphaDose{1}(V,:)*wOnes)) / (sum((dij.mSqrtBetaDose{1}(V,:) * wOnes).^2));
           q            = -(effectTarget * length(V)) / (sum((dij.mSqrtBetaDose{1}(V,:) * wOnes).^2));
           wInit        = -(p/2) + sqrt((p^2)/4 -q) * wOnes;

    elseif isequal(pln.bioOptimization,'LSM_RBExD') ||isequal(pln.bioOptimization,'LEMIV_RBExD')
        
           %pre-calculations
           dij.gamma      = zeros(dij.numOfVoxels,1);
           idx            = dij.bx~=0; 
           dij.gamma(idx) = dij.ax(idx)./(2*dij.bx(idx)); 
            
           % calculate current in target
           CurrEffectTarget = (dij.mAlphaDose{1}(V,:)*wOnes + (dij.mSqrtBetaDose{1}(V,:)*wOnes).^2);
           % ensure a underestimated biological effective dose 
           TolEstBio        = 1.2;
           % calculate maximal RBE in target
           maxCurrRBE = max(-cst{ixTarget,5}.alphaX + sqrt(cst{ixTarget,5}.alphaX^2 + ...
                        4*cst{ixTarget,5}.betaX.*CurrEffectTarget)./(2*cst{ixTarget,5}.betaX*(dij.physicalDose{1}(V,:)*wOnes)));
           wInit      =  ((doseTarget)/(TolEstBio*maxCurrRBE*max(dij.physicalDose{1}(V,:)*wOnes)))* wOnes;
    end
    
else 
    bixelWeight         =  (doseTarget)/(mean(dij.physicalDose{1}(V,:)*wOnes)); 
    wInit               = wOnes * bixelWeight;
    pln.bioOptimization = 'none';
end

%% ToDo: incorporate changes of robOpt in GUI

% check if deterministic / stoachastic optimization is turned on
if isfield(pln,'robOpt')
   if strcmp(pln.robOpt,'none')
       dij.indexforOpt    = 1;
       dij.numOfScenarios = 1;
   else
       dij.indexforOpt    = find(~cellfun(@isempty, dij.physicalDose))'; 
       dij.numOfScenarios = numel(dij.indexforOpt);
   end
else
      pln.robOpt = 'none';
      dij.indexforOpt    = 1;
      dij.numOfScenarios = 1;
end

% set optimization options
options.radMod          = pln.radiationMode;
options.bioOpt          = pln.bioOptimization;
options.robOpt          = pln.robOpt;
options.ID              = [pln.radiationMode '_' pln.bioOptimization];
options.numOfScenarios  = dij.numOfScenarios;

% set callback functions.
funcs.objective         = @(x) matRad_objFuncWrapper(x,dij,cst,options);
funcs.constraints       = @(x) matRad_constFuncWrapper(x,dij,cst,options);
funcs.gradient          = @(x) matRad_gradFuncWrapper(x,dij,cst,options);
funcs.jacobian          = @(x) matRad_jacobFuncWrapper(x,dij,cst,options);
funcs.jacobianstructure = @( ) matRad_getJacobStruct(dij,cst);

% scale objective and constraint function
gInit    = abs(matRad_gradFuncWrapper(wInit,dij,cst,options));
fScaling = 1e2/max(gInit);


if ~isempty(matRad_getConstBoundsWrapper(cst,options))

    jInit    = abs(matRad_jacobFuncWrapper(wInit,dij,cst,options));
        
    for i = 1:length(matRad_getConstBoundsWrapper(cst,options))  
        
        wInitTmp = wInit;
        while sum(jInit(i,:)) == 0
            wInitTmp = wInitTmp - 0.1*wInit;
            jInit = abs(matRad_jacobFuncWrapper(wInitTmp,dij,cst,options));
        end
        cScalingTmp(i,1) = 1e-3/max(jInit(i,:));
    end
    cScaling = cScalingTmp;
end

options.ipopt.acceptable_constr_viol_tol = max(cScaling)*options.ipopt.acceptable_constr_viol_tol;
[options.cl,options.cu] = matRad_getConstBoundsWrapper(cst,options);  

% Run IPOPT.
[wOpt, info] = ipopt(wInit,funcs,options);

% calc dose and reshape from 1D vector to 2D array
fprintf('Calculating final cubes...\n');

resultGUI = matRad_calcCubes(wOpt,dij,cst,1);
resultGUI.wUnsequenced = wOpt;

% save optimization info in resultGUI
resultGUI.optInfo.IPOPTinfo                                 = info;
resultGUI.optInfo.IPOPToptions                              = options;
resultGUI.optInfo.wInit                                     = wInit;
resultGUI.optInfo.wOpt                                      = wOpt;
resultGUI.optInfo.globalVar.matRad_global_x                 = matRad_global_x;
resultGUI.optInfo.globalVar.matRad_global_d                 = matRad_global_d;
resultGUI.optInfo.globalVar.matRad_STRG_C_Pressed           = matRad_STRG_C_Pressed;
resultGUI.optInfo.globalVar.matRad_objective_function_value = matRad_objective_function_value;
resultGUI.optInfo.globalVar.matRad_iteration                = matRad_iteration;
resultGUI.optInfo.globalVar.kDVH                            = kDVH;
resultGUI.optInfo.globalVar.kDCH                            = kDCH;
resultGUI.optInfo.globalVar.JACOBIAN                        = JACOBIAN;
resultGUI.optInfo.globalVar.GRADIENT                        = GRADIENT;
resultGUI.optInfo.globalVar.CONSTRAINT                      = CONSTRAINT;
resultGUI.optInfo.globalVar.fScaling                        = fScaling;
resultGUI.optInfo.globalVar.cScaling                        = cScaling;

% unset Key Pressed Callback of Matlab command window
if ~isdeployed
    set(h_cw, 'KeyPressedCallback',' ');
end

% clear global variables
clearvars -global matRad_global_x matRad_global_d matRad_objective_function_value matRad_STRG_C_Pressed matRad_iteration;
clearvars -global kDVH kDCH GRADIENT JACOBIAN fScaling cScaling CONSTRAINT

% unblock mex files
clear mex