function [D, pln, dij, resultGUI] = matRad_4dRBEcalc(ct, cst, File)

% calculates all variations of 4D and RBE calculations
%  
%
% input
%   ct:          matRads ct struct
%   cst:         matRads cst struct
%   pln:         matRads pln struct
%   dij:         inclusive alpha beta cubes
%   resultGUI    
%
%   
addpath('D:\Matrad\')
addpath('D:\Matrad\4Ddose')

% output D with all dose cubes
%D.name = {'Dopt', 'Drecalc3D', 'Dopt -Drecalc3D', 'Drecalc4Dconst', 'Drecalc4Dvar', 'Drecalc4Dconst-Drecalc4Dvar','Dopt - Drecalc4Dconst', 'Drecalc3D-Drecalc4Dvar', 'Dopt-Drecalc4Dvar'};
D.name = {'(A)', '(B)', '(A)-(B)', '(C)', '(D)', '(C)-(D)','(A)-(C)', '(B)-(D)', '(A)-(D)'};

D.isolines = {1, 1, 0, 1, 1, 0, 0, 0, 0};


%Pln info
pln.numOfFractions  = 30;
D.fractions = pln.numOfFractions;  
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'HITfixedBL';

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longSpotSpacing = 5;      % only relevant for HIT machine, not generic
pln.propStf.gantryAngles    = [210]; % [?] ;
pln.propStf.couchAngles     = [0 ]; % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%optimization settings
pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
modelName    = 'MCN';             % none: for photons, protons, carbon            % constRBE: constant RBE 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          


% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);

%calc dose matrix
stf = matRad_generateStf(ct,cst,pln);

dij = matRad_calcParticleDose(ct,stf,pln,cst,false);

%opt for const RBE
modelName = 'constRBE';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf);
Dopt = resultGUI.RBExD;
matRad_export_HITXMLPlan_modified('PlanRBE_new',  pln, stf, resultGUI, 'stfMode')  

%% break for makeLmdout

%4D for const RBE
[resultGUI, ~] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI,  File);

Drecalc4Dconst = resultGUI.accRBExD;

%4D for variable RBE
modelName = 'MCN';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
[resultGUI, ~] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI,  File);
Drecalc4Dvar = resultGUI.accRBExD;

%recalc for variable RBE
Drecalc3D  = matRad_calcMcNRBExD(dij, cst, resultGUI);

D.data = {Dopt, Drecalc3D, Dopt-Drecalc3D, Drecalc4Dconst, Drecalc4Dvar, Drecalc4Dconst-Drecalc4Dvar, Dopt - Drecalc4Dconst, Drecalc3D-Drecalc4Dvar, Dopt-Drecalc4Dvar};






