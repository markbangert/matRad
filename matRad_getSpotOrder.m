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

if nargin < 2
    plotting = 'off';
end

order = zeros(sum([stf.totalNumOfBixels]), 1);
bixelInfo = struct;

% first loop loops over all bixels to store their position and ray number
% in each IES
wOffset = 0;
for i = 1:length(stf) % looping over all beams
    
    usedEnergies = unique([stf(i).ray(:).energy]);
    usedEnergiesSorted = sort(usedEnergies, 'descend');
        
    for e = 1:length(usedEnergies) % looping over IES's
        
        s = 1;
        
        for j = 1:stf(i).numOfRays % looping over all rays
            
            % find the rays which are active in current IES
            if(any(stf(i).ray(j).energy == usedEnergiesSorted(e)))
                
                x = stf(i).ray(j).rayPos_bev(1);
                y = stf(i).ray(j).rayPos_bev(3);
                
                bixelInfo(i).IES(e).x(s)       = x; % store x position
                bixelInfo(i).IES(e).y(s)       = y; % store y position
                bixelInfo(i).IES(e).w_index(s) = wOffset + ...
                                                 sum(stf(i).numOfBixelsPerRay(1:(j-1))) + ...
                                                 find(stf(i).ray(j).energy == usedEnergiesSorted(e)); % store index
                
                s = s + 1;
                
            end
        end
    end
    
    wOffset = wOffset + sum(stf(i).numOfBixelsPerRay);
    
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
        x_sorted = sort(bixelInfo(i).IES(e).x, 'ascend');
        
        for k = 1:length(y_sorted)
            
            y = y_sorted(k);
            x = x_sorted(k);
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
                
                x = bixelInfo(i).IES(e).x(ss);
                
                w_index = bixelInfo(i).IES(e).w_index(ss);
                
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
                order_counter  = order_counter + 1;
            end
        end
    end
end

end
