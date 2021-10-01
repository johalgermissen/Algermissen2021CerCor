function output = taft_prune(input)

% output = taft_prune(input)
% 
% Prune away clusters that have "bridges" of only 1 voxel.
% 
% INPUTS:
% input     = 2D matrix with 1 for cluster above threshold, otherwise 0.
% OUTPUTS:
% output    = same matrix, but bridges pruned, so now more smaller
% clusters.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

%% Initialize:

output = input; % initialize output object

%% Prune isolated voxels (i.e. bridges) by evaluating connectivity

for xPos = 1:size(output,1) % for each column in 2D matrix
    for yPos = 1:size(output,2) % for each row in 2D matrix
        
        % Only for voxels with value of 1:
        % Only for voxels not at the edges of the matrix:
        if (output(xPos,yPos) == 1) && (xPos > 1 && xPos < size(output,1) && yPos > 1 && yPos < size(output,2))
            
            %% 1) Voxels with only two connections (i.e. diagonal edges):
            
            iEdge = 0;

            %% A) Alternative A: loop over X and Y coordinates, drop middle of grid:
            
            % Loop through all points in 3x3 grid with voxel itself in the middle:
            for xxPos = [xPos-1,xPos,xPos+1] % x positions
                for yyPos = [yPos-1,yPos,yPos+1] % y positions
                    if ~(xxPos == xPos && yPos == yPos) % if not middle of grid itself, but edge of grid:
                        if output(xxPos,yyPos) == 1 % if edge itself 1:
                            iEdge = iEdge + 1; % count edges
                        end % end if output is 1
                    end % end if not middle
                end % end yyPos
            end % end xxPos
            
            %% B) Alternative B: consider each possible edge separately manually:
            
%             % Left connection:
%             if output(xPos-1,yPos) == 1
%                 iEdge = iEdge + 1;
%             end
%             % Right connection:
%             if output(xPos+1,yPos) == 1
%                 iEdge = iEdge + 1;
%             end
%             % Above connection:
%             if output(xPos,yPos-1) == 1
%                 iEdge = iEdge + 1;
%             end
%             % Below connection:
%             if output(xPos,yPos+1) == 1
%                 iEdge = iEdge + 1;
%             end
%             % Left above connection:
%             if output(xPos-1,yPos-1) == 1
%                 iEdge = iEdge + 1;
%             end
%             % Left below connection:
%             if output(xPos+1,yPos-1) == 1
%                 iEdge = iEdge + 1;
%             end
%             % Right above connection:
%             if output(xPos-1,yPos+1) == 1
%                 iEdge = iEdge + 1;
%             end
%             % Right below connection:
%             if output(xPos+1,yPos+1) == 1
%                 iEdge = iEdge + 1;
%             end
            
            %% C) Evaluate if only 2 edges or less:
            
            if iEdge <= 2 % if only two neighbors (i.e. diagonal bridge)
                output(xPos,yPos) = 0; % remove bridge
            end
            
            %% 2) Prune only horizontal bridges:
            
%             if output(xPos-1,yPos) == 1 && output(xPos+1,yPos) == 1 && output(xPos,yPos-1) == 0 && output(xPos,yPos+1) == 0
%                 output(xPos,yPos) = 0; % remove bridge
%             end
            
            %% 3) Prune only vertical bridges:
            
%             if output(xPos-1,yPos) == 0 && output(xPos+1,yPos) == 0 && output(xPos,yPos-1) == 1 && output(xPos,yPos+1) == 1
%                 output(xPos,yPos) = 0; % remove bridge
%             end
            
        end % end if middle is 1 and not edge of entire matrix
    end % end iCol
end % end iRow 

end % end of function.