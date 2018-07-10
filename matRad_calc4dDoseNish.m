function [resultGUI, bixelInfo, w] = matRad_calc4dDoseNish(ct, pln, dij, stf, cst, resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad 4D dose calculation
% 
% call
%   
%
% input 
%   
%
% output
%   
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a time sequence for when each bixel is irradiated, the sequence
% follows the backforth spot scanning
bixelInfo = matRad_makeBixelTimeSeq(stf, resultGUI.w);

motionPeriod = 5; % a whole breathing motion period (in seconds)
motion = 'linear'; % the assumed motion of chest 
numOfPhases = size(ct.cube, 2);

% prepare a phase matrix in which place each bixel dose in it's phase
bixelInfo = matRad_makePhaseMatrix(bixelInfo, numOfPhases, motionPeriod, motion);
      