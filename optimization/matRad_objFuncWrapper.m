function f = matRad_objFuncWrapper(w,dij,cst,type)

global matRad_voxelWeighting;

% initialize voxel calc Flag
[matRad_voxelWeighting{:,2}] = deal(true);

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

% Initialize f
f = 0;

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})

            % only perform gradient computations for objectives
            if isempty(strfind(cst{i,6}(j).type,'constraint'))

                % compute reference
                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) &&...
                    isequal(type,'effect') 

                    d_ref = dij.ax(cst{i,4}{1}).*cst{i,6}(j).dose + dij.bx(cst{i,4}{1})*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
               % if conventional opt: just sum objectives of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref);

                % if prob opt: sum up expectation value of objectives
                elseif strcmp(cst{i,6}(j).robustness,'probabilistic')

                    for k = 1:dij.numOfScenarios

                        d_i = d{k}(cst{i,4}{1});

                        f = f + dij.probOfScenarios(k) * matRad_objFunc(d_i,cst{i,6}(j),d_ref);

                    end

                % if voxel-wise worst case: sum up objective of min/max dose
                elseif strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')

                    % prepare min/max dose vector we have chosen voxel-wise worst case
                    if ~exist('d_max','var')
                        d_max = max([d{:}],[],2);
                        d_min = min([d{:}],[],2);
                    end

                    if isequal(cst{i,3},'OAR')
                        d_i = d_max(cst{i,4}{1});
                    elseif isequal(cst{i,3},'TARGET')
                        d_i = d_min(cst{i,4}{1});
                    end

                    f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                    
                % if coveraged based opt    
                elseif strcmp(cst{i,6}(j).robustness,'coverage')
                    
                    d_i = [];
                    
                    % get cst index of VOI that corresponds to VOI ring
                    cstidx = find(strcmp(cst(:,2),cst{i,2}(1:end-4)));
                    
                    % get dose of VOI that corresponds to VOI ring
                    for k = 1:dij.numOfScenarios
                        d_i{k} = d{k}(cst{cstidx,4}{1});
                    end
                    
                    % calc invers DCH of VOI
                    refQ   = cst{i,6}(j).coverage/100;
                    refVol = cst{i,6}(j).volume/100;
                    d_ref2 = matRad_calcInversDCH(refVol,refQ,d_i,dij.numOfScenarios);
                    
                    % get dose of VOI ring
                    d_i = d{1}(cst{i,4}{1});
                    
                    % get voxel dependetn weigthing
                    if isequal(cst{i,5}.voxelWeightingType,'heurWeighting')
                        matRad_calcVoxelWeighting(i,j,cst,d_i,d_ref,d_ref2)
                        voxelWeighting = matRad_voxelWeighting{i,1};
                        
                    elseif isequal(cst{i,5}.voxelWeightingType,'probWeighting')
                        voxelWeighting = 5*cst{i,5}.voxelProb;    
                    end
    
                    f = f + matRad_objFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);

                end
            
            end
       
        end
            
    end
    
end
