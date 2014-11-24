function P = ComputeTransitionProbabilitiesI( stateSpace, controlSpace, disturbanceSpace, mazeSize, walls, targetCell )
%COMPUTETRANSITIONPROBABILITIESI Compute transition probabilities.
% 	Compute the transition probabilities between all states in the state
%   space for all attainable control inputs.
%
%   P = ComputeTransitionProbabilitiesI(stateSpace, controlSpace,
%   disturbanceSpace, mazeSize, walls, targetCell) computes the transition
%   probabilities between all states in the state space for all attainable
%   control inputs.
%
%   Input arguments:
%
%       stateSpace:
%           A (MN x 2) matrix, where the i-th row represents the i-th
%           element of the state space. Note that the state space also
%           contains the target cell, in order to simplify state indexing.
%
%       controlSpace:
%           A (L x 2) matrix, where the l-th row represents the l-th
%           element of the control space.
%
%       disturbanceSpace:
%           A (S x 3) matrix 'disturbanceSpace', where the first two
%           columns of each row represent an element of the disturbance
%           space, and the third column represents the corresponding
%           probability.
%
%       mazeSize:
%           A (1 x 2) matrix containing the width and the height of the
%           maze in number of cells.
%
%   	walls:
%          	A (2 x 2K) matrix containing the K wall segments, where the start
%        	and end point of the k-th segment are stored in column 2k-1
%         	and 2k, respectively.
%
%    	targetCell:
%          	A (2 x 1) matrix describing the position of the target cell in
%         	the maze.
%
%   Output arguments:
%
%       P:
%           A (MN x MN x L) matrix containing the transition probabilities
%           between all states in the state space for all attainable
%           control inputs. The entry P(i, j, l) represents the transition
%           probability from state i to state j if control input l is
%           applied.

MN = size(stateSpace,1);
M = mazeSize(1);
L = size(controlSpace,1);
P = zeros(MN,MN,L);
wallsMatrix = GenerateWallsMatrix(mazeSize, walls);
for cell = 1:MN
    controls = availableControls(cell,wallsMatrix);
    for i = 1:size(controls)
        [x,y] = controlSpace(controls(i));
        cellArrival = cell + x*M + y;
        wallsCellArrival = wallsMatrix(cellArrival);
        
        P(cell,cellArrival,controls(i)) = 1/5;
        %RIGHT
        if(wallsCellArrival(1) == 0)
            P(cell,cellArrival + M,controls(i)) = 1/5;
        else
            P(cell,cellArrival,controls(i)) = P(cell,cellArrival,controls(i)) + 1/5;
        end
        %UP
        if(wallsCellArrival(2) == 0)
            P(cell,cellArrival + 1,controls(i)) = 1/5;
        else
            P(cell,cellArrival,controls(i)) = P(cell,cellArrival,controls(i)) + 1/5;
        end
        %LEFT
        if(wallsCellArrival(3) == 0)
            P(cell,cellArrival - M,controls(i)) = 1/5;
        else
            P(cell,cellArrival,controls(i)) = P(cell,cellArrival,controls(i)) + 1/5;
        end
        %BOTTOM
        if(wallsCellArrival(4) == 0)
            P(cell,cellArrival - 1,controls(i)) = 1/5;
        else
            P(cell,cellArrival,controls(i)) = P(cell,cellArrival,controls(i)) + 1/5;
        end 
    end
end
end

%% FUNCTION TO FIND AVAILABLE CONTROLS FOR ONE SPECIFIC CELL
function controls = availableControls(cell,wallsMatrix)
controls = 7;
%RIGHT CONTROLS
if (wallsMatrix(cell,1) == 0)
   controls = [controls, 11];
   cellRight = cell + m; 
   if(wallsMatrix(cellRight,1) == 0)
      controls = [controls, 13]; 
   end
end
%UP CONTROLS
if (wallsMatrix(cell,2) == 0)
   controls = [controls, 8];
   cellUp = cell + 1; 
   if(wallsMatrix(cellUp,2) == 0)
      controls = [controls, 9]; 
   end
end
%LEFT CONTROLS
if (wallsMatrix(cell,3) == 0)
   controls = [controls, 3];
   cellLeft = cell + m; 
   if(wallsMatrix(cellLeft,3) == 0)
      controls = [controls, 1]; 
   end
end
%BOTTOM CONTROLS
if (wallsMatrix(cell,4) == 0)
   controls = [controls, 6];
   cellBottom = cell - 1; 
   if(wallsMatrix(cellBottom,4) == 0)
      controls = [controls, 5]; 
   end
end
%DIAGONALS CONTROLS RIGHT-UP
if (wallsMatrix(cell,1) == 0 && wallsMatrix(cell,2) == 0)
   cellRight = cell + m;
   cellUp = cell + 1;
   if(wallsMatrix(cellRight,2) == 0 && wallsMatrix(cellUp,1) == 0)
      controls = [controls, 12]; 
   end
end
%DIAGONALS CONTROLS UP-LEFT
if (wallsMatrix(cell,2) == 0 && wallsMatrix(cell,3) == 0)
   cellUp = cell + 1;
   cellLeft = cell - m;
   if(wallsMatrix(cellUp,3) == 0 && wallsMatrix(cellLeft,2) == 0)
      controls = [controls, 4]; 
   end
end
%DIAGONALS CONTROLS LEFT-BOTTOM
if (wallsMatrix(cell,3) == 0 && wallsMatrix(cell,4) == 0)
   cellLeft = cell - m;
   cellBottom = cell - 1;
   if(wallsMatrix(cellLeft,4) == 0 && wallsMatrix(cellBottom,3) == 0)
      controls = [controls, 2]; 
   end
end
%DIAGONALS CONTROLS BOTTOM-RIGHT
if (wallsMatrix(cell,4) == 0 && wallsMatrix(cell,1) == 0)
   cellBottom = cell - 1;
   cellRight = cell + m;
   if(wallsMatrix(cellBottom,1) == 0 && wallsMatrix(cellRight,4) == 0)
      controls = [controls, 10]; 
   end
end
end
