function cst = matRad_coverageBasedCstManipulation(cst,ct,multScen,normalTissueRingVOIName,varargin)

covFlag = 0;
Counter = 0;

for  i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        if sum(strcmp({cst{i,6}(:).robustness},'coverage')) > 0 & sum(strcmp({cst{i,6}(:).type},'min DCH objective')) > 0 |...
           sum(strcmp({cst{i,6}(:).robustness},'coverage')) > 0 & sum(strcmp({cst{i,6}(:).type},'max DCH objective')) > 0
       
           if isempty(strfind(cst{i,2},'Union'))
                
                covFlag = 1;
                Counter = Counter + 1;  

                % create VOI scenario union structure
                cstVOIScenUnion{Counter,1} = size(cst,1) - 1 + Counter;
                cstVOIScenUnion{Counter,2} = [cst{i,2},' ScenUnion'];
                cstVOIScenUnion{Counter,3} = cst{i,3};
                cstVOIScenUnion{Counter,5} = cst{i,5};
                
                
                if multScen.numOfShiftScen > 1
                    % use dij scenarios

                    [yCoordsVOI_vox, xCoordsVOI_vox, zCoordsVOI_vox] = ind2sub(ct.cubeDim,cst{i,4}{1});
                    voxelProbCube                                    = zeros(ct.cubeDim);
                    VOIScenUnionidx                                  = cst{i,4}{1};

                    for k = 1:multScen.numOfShiftScen

                        % round shifts to voxel dimensions
                        xShift_vox = round(multScen.shifts(1,k)/ct.resolution.x);
                        yShift_vox = round(multScen.shifts(2,k)/ct.resolution.y);
                        zShift_vox = round(multScen.shifts(3,k)/ct.resolution.z);

                        shiftedVOIidx = sub2ind(ct.cubeDim, yCoordsVOI_vox - yShift_vox,...
                                                            xCoordsVOI_vox - xShift_vox,...
                                                            zCoordsVOI_vox - zShift_vox);

                        % fill prob cube                                   
                        voxelProbCube(shiftedVOIidx) = voxelProbCube(shiftedVOIidx) + 1/multScen.numOfShiftScen;                                   

                        % build scenario union structure
                        VOIScenUnionidx = union(VOIScenUnionidx,shiftedVOIidx);

                    end

                    cstVOIScenUnion{Counter,4}{1}        = VOIScenUnionidx;
                    cstVOIScenUnion{Counter,5}.voxelProb = voxelProbCube(cstVOIScenUnion{Counter,4}{1})';

                elseif multScen.numOfShiftScen == 1
                    % create scnearios with shifts

                    % sample VOI shifts
                    cstVOIScenUnion{Counter,5}.VOIShift = matRad_sampleVOIShift(cst,ct,multScen.shiftSize,cst{i,2},1000);
                    cstVOIScenUnion{Counter,4}{1}       = find(cstVOIScenUnion{Counter,5}.VOIShift.voxelProbCube > 0);

                    for k  = 1:size(cst,1)
                        cst{k,5}.VOIShift = cstVOIScenUnion{Counter,5}.VOIShift;
                    end  

                    cstVOIScenUnion{Counter,5}.VOIShift.voxelProb = cstVOIScenUnion{Counter,5}.VOIShift.voxelProbCube(cstVOIScenUnion{Counter,4}{1});

                end                                    

                % pass coverage based objective specification to VOI scenario union structure
                logidx                     = strcmp({cst{i,6}(:).robustness},'coverage');
                cstVOIScenUnion{Counter,6} = cst{i,6}(logidx);
                if sum(~logidx) == 0
                    cst{i,6} = []; 
                else
                    cst{i,6} = cst{i,6}(~logidx);    
                end
               
           end
        end 
    end
end

if covFlag
    cst = [cst;cstVOIScenUnion];
end

if ischar(normalTissueRingVOIName)  
    cstidx = find(strcmp([cst(:,2)],'BODY'));
    
    cstTmp{1,1} = cst{end,1} + 1;
    cstTmp{1,2} = 'normal tissue Ring';
    cstTmp{1,3} = cst{cstidx,3};
    cstTmp{1,5} = cst{cstidx,5};
    cstTmp{1,6} = cst{cstidx,6};
    cst{cstidx,6} = [];
    
    voiCube                                                             = zeros(ct.cubeDim);
    voiCube(cst{find(strcmp([cst(:,2)],normalTissueRingVOIName)),4}{1}) = 1;
    
    % inner margin
    [margin.x, margin.y, margin.z] = deal(varargin{1});
    voiCubeInner = matRad_addMargin(voiCube,cst,ct.resolution,margin,true);
    VInner    = find(voiCubeInner>0);

    % outer margin
    [margin.x, margin.y, margin.z] = deal(varargin{2});
    voiCubeOuter = matRad_addMargin(voiCube,cst,ct.resolution,margin,true);
    VOuter    = find(voiCubeOuter>0);
    
    cstTmp{1,4}{1} = setdiff(VOuter,VInner);
    
    cst = [cst;cstTmp];
end

end