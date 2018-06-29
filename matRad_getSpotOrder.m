function order = matRad_getSpotOrder(stf, plotting)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad inverse planning wrapper function
%
% call
%   [resultGUI,info] = matRad_fluenceOptimization(dij,cst,pln)
%
% input
%   stf:            matRad steering information struct
%   plotting:       sets the plotting option 'on' or 'off' (temporary)
%
% output
%   order:          order of each bixel according to the spot scanning
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order = zeros(sum([stf.totalNumOfBixels]), 1);
bixelInfo = struct;


% first loop loops over all bixels to store their position and ray number
% in each IES
for i = 1:length(stf) % looping over all beams
    
    usedEnergies = unique([stf(i).ray(:).energy]);
    usedEnergiesSorted = sort(usedEnergies, 'descend');
    
    
    for e = 1: length(usedEnergies) % looping over IES's
        s = 1;
        for j = 1:stf(i).numOfRays % looping over all rays
            
            % find the rays which are active in current IES
            if(any(stf(i).ray(j).energy == usedEnergiesSorted(e)))
                
                x = stf(i).ray(j).rayPos_bev(1);
                y = stf(i).ray(j).rayPos_bev(3);
                
                bixelInfo(i).IES(e).x(s) = x; % store x position
                bixelInfo(i).IES(e).y(s) = y; % store y position
                bixelInfo(i).IES(e).j(s) = j; % store ray number
                
                s = s + 1;
                
            end
        end
    end
end

% after storing all the required information,
% same loop over all bixels will put each bixel in it's order

order_counter = 1;

for i = 1:length(stf)
    
    usedEnergies = unique([stf(i).ray(:).energy]);
    usedEnergiesSorted = sort(usedEnergies, 'descend');
    temp = 0; % temporary variable for plotting figures
    
    for e = 1: length(usedEnergies)
        
        % sort the y positions from high to low (backforth is up do down)
        y_sorted = sort(unique(bixelInfo(i).IES(e).y), 'descend');
        
        for k = 1:length(y_sorted)
            
            y = y_sorted(k);
            % find indexes corresponding to current y position
            % in other words, number of bixels in the current row
            index = find(bixelInfo(i).IES(e).y == y);
            
            % since backforth fasion is zig zag like, flip the order every
            % second row
            if ~rem(k,2)
                index = fliplr(index);
            end
            
            % loop over all the bixels in the row
            for ss = index
                
                j = bixelInfo(i).IES(e).j(ss);
                x = bixelInfo(i).IES(e).x(ss);
                
                % follwing if block assigns corresponding stf index for the
                % current ray dependent on (beam number, ray number, and
                % energy slice)
                if (i ~= 1) && (j ~= 1) %to avoid zero indexing
                    w_index = (i - 1) * stf(i-1).totalNumOfBixels ...
                        + sum(stf(i).numOfBixelsPerRay(1:j-1)) ...
                        + find(stf(i).ray(j).energy == usedEnergiesSorted(e));
                    
                    % this is how w_index is calculated:
                    % w_index = number of bixels in previous beam
                    %         + number of bixels in previous rays
                    %         + number of previous bixels in current ray
                    
                    % following does the same as above, they are separated to
                    % avoid zero indexing
                elseif (j == 1) && (i ~= 1)
                    w_index = (i - 1) * stf(i-1).totalNumOfBixels ...
                        + find(stf(i).ray(j).energy == usedEnergiesSorted(e));
                elseif (i == 1) && (j~=1)
                    w_index = sum(stf(i).numOfBixelsPerRay(1:j-1)) ...
                        + find(stf(i).ray(j).energy == usedEnergiesSorted(e));
                else
                    w_index = find(stf(i).ray(j).energy == usedEnergiesSorted(e));
                end
                
                % following is to plot the bixel ordering, can be deleted
                % in the final version of the code
                if(strcmp(plotting, 'on'))
                    
                    if(temp ~= e)
                        clf
                        h = animatedline('LineStyle', 'none', 'Marker', 'o');
                        axis([-50 50 -50 50])
                        title(['Beam #', num2str(i), ', IES #', num2str(e)])
                    end
                    
                    temp = e;
                    
                    addpoints(h, x, y);
                    strmin = [num2str(order_counter), '  '];
                    text(x,y,strmin,'HorizontalAlignment','right');
                    pause(.1);
                    drawnow
                end
                
                % assign the order to the corresponding stf index
                order(w_index) = order_counter;
                order_counter = order_counter + 1;
            end
        end
    end
end

end
