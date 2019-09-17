function [ctSample, cstSample] = lstm_sampleGeo(origCt,origCst,isoCenter)

resolution = 2; % [mm]
lateralDim = 11;
longitudinalDim = 300/resolution;
lateralSamplingPoints = resolution * [-lateralDim:lateralDim];
longitudinalSamplingPoints = resolution * [-longitudinalDim:longitudinalDim];

% set up a meshgrid for extraction of sampling points
[X,Y,Z] = meshgrid(lateralSamplingPoints,longitudinalSamplingPoints,lateralSamplingPoints);

reshapeDim = size(X);

% sample superior gantry angles
ga = mod(270+rand*180,360);

% sample offset in z direction
zOffset = (2*rand-1) * 30;

% rotate extraction points and add offset
R = matRad_getRotationMatrix(ga,0);
rotCoords = [X(:) Y(:) Z(:)] * R;
rotCoords(:,1) = rotCoords(:,1) + isoCenter(1);
rotCoords(:,2) = rotCoords(:,2) + isoCenter(2);
rotCoords(:,3) = rotCoords(:,3) + isoCenter(3) + zOffset;

ct = interp3(origCt.resolution.x*[1:origCt.cubeDim(1)], ...
                   origCt.resolution.y*[1:origCt.cubeDim(2)], ...
                   origCt.resolution.z*[1:origCt.cubeDim(3)],origCt.cubeHU{1}, ...
                   rotCoords(:,1), ...
                   rotCoords(:,2), ...
                   rotCoords(:,3), ...
                   'cubic');
               
ct = reshape(ct,reshapeDim);

% find first slice with value larger than threshold
thresholdHU = -750;
for i = 1:size(ct,1)
    if any(ct(i,:) > thresholdHU)
        break
    end
end

% pad or crop ct so that it always has the same dimension...
lengthOfCt_mm = 300; % [mm]
lengthOfCt_vox = lengthOfCt_mm/resolution;

% cut start
ct = ct(i:end,:,:);

% crop or pad end
if size(ct,1) > lengthOfCt_vox
    ct = ct(1:lengthOfCt_vox,:,:);
elseif size(ctSamplel,1) < lengthOfCt_vox
    ct(lengthOfCt_vox,size(ct,2),size(ct,3)) = NaN;
end

% replace NaNs
ct(isnan(ct)) = -1024;

% % for visual debugging
%imagesc(squeeze(ct(:,:,8)))
%axis equal tight
%drawnow

% prepare matRad ct structure
ctSample.cubeDim = size(ct);
ctSample.cubeHU{1} = ct;
ctSample.resolution.x = resolution;
ctSample.resolution.y = resolution;
ctSample.resolution.z = resolution;
ctSample.numOfCtScen = 1;
ctSample.x = resolution * [1:ctSample.cubeDim(2)];
ctSample.y = resolution * [1:ctSample.cubeDim(1)];
ctSample.z = resolution * [1:ctSample.cubeDim(3)];

% dummy Cst
cstSample = origCst(1,1:6);
cstSample{1,2} = 'dummy';
cstSample{1,3} = 'TARGET';
cstSample{1,4}{1} = 1:numel(ct); 
