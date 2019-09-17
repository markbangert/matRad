% matRad script
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%matRad_rc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
%load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat

% meta information for treatment plan

pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'generic_MCsquare';

pln.numOfFractions  = 1;

% beam geometry settings
pln.propStf.bixelWidth      = 500; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [0]; % [?]
pln.propStf.couchAngles     = [0]; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 2; % [mm]

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

% generate steering file
stf = matRad_generateStf(ct,cst,pln);

stf.totalNumOfBixels = 1;
stf.numOfBixelsPerRay = 1;
stf.ray.energy(1:end-1) = [];
stf.ray.focusIx(1:end-1) = [];
stf.ray.rangeShifter(1:end-1) = [];
stf.ray.weight = 1;

%% sample dose calculation
% sample a geometry
[ctSample, cstSample] = lstm_sampleGeo(ct,cst,pln.propStf.isoCenter);

% adjust isocenter
stf.isoCenter = [ctSample.resolution.x ctSample.resolution.y ct.resolution.z] .* (ctSample.cubeDim([2 1 3])-1)/2;

% calculate analytical and MC dose
doseMC = matRad_calcParticleDoseMC(ctSample,stf,pln,cstSample,100000,1);
doseAna = matRad_calcParticleDose(ctSample,stf,pln,cstSample,1);

% reshape result to cube
doseCubeMC = reshape(full(doseMC.physicalDose{1}),ctSample.cubeDim);
doseCubeAna = reshape(full(doseAna.physicalDose{1}),ctSample.cubeDim);

% visualization
slice = (ctSample.cubeDim(3)-1)/2;

subplot(1,4,1)
imagesc(ctSample.cubeHU{1}(:,:,slice))
axis equal tight
colorbar
title('CT')

subplot(1,4,2)
imagesc(doseCubeMC(:,:,slice))
axis equal tight
colorbar
title('MC dose')

subplot(1,4,3)
imagesc(doseCubeAna(:,:,slice))
axis equal tight
colorbar
title('Analytical dose')

subplot(1,4,4)
imagesc(doseCubeMC(:,:,slice)-doseCubeAna(:,:,slice))
axis equal tight
colorbar
title('MC - Analytical dose')
