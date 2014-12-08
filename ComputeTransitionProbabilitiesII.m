function P = ComputeTransitionProbabilitiesII( stateSpace, controlSpace, disturbanceSpace, mazeSize, walls, targetCell, holes, resetCell )
%COMPUTETRANSITIONPROBABILITIESII Compute transition probabilities.
% 	Compute the transition probabilities between all states in the state
%   space for all attainable control inputs.
%
%   P = ComputeTransitionProbabilitiesII(stateSpace, controlSpace,
%   disturbanceSpace, mazeSize, walls, targetCell, holes, resetCell)
%   computes the transition probabilities between all states in the state
%   space for all attainable control inputs.
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
%    	holes:
%         	A (2 x H) matrix containg the H holes of the maze. Each column
%         	represents the position of a hole.
%
%   	resetCell:
%         	A (2 x 1) matrix describing the position of the reset cell in
%           the maze.
%
%   Output arguments:
%
%       P:
%           A (MN x MN x L) matrix containing the transition probabilities
%           between all states in the state space for all attainable
%           control inputs. The entry P(i, j, l) represents the transition
%           probability from state i to state j if control input l is
%           applied.

%inizialize some useful dimension
MN = size(stateSpace,1);
M = mazeSize(2);
L = size(controlSpace,1);
P = zeros(MN,MN,L);

resetCell = (resetCell(1) - 1)*M + resetCell(2);
%create the matrix of the HOLES
holeSpace = zeros(1,size(holes,2));
for i = 1 : size(holes,2)
   holeSpace(i) = (holes(1,i) - 1)*M + holes(2,i);
end

%create the matrix of the WALLS
wallsMatrix = GenerateWallsMatrix(mazeSize, walls);

for cell = 1:MN
    if(cell == (targetCell(1)-1)*M +targetCell(2)) %if cell is target
        %only control aviable:STAY
        P(cell,cell,7) = 1;
    else
        %check all APPLICABLE CONTROLS
        controls = applicableControls(cell,wallsMatrix,M,holeSpace);
        for i = 1:size(controls,2)
            coords = controlSpace(controls(i),:);
            %apply CONTROL
            cellArrWithCont = cell + coords(1)*M + coords(2);
            if(ismember(cellArrWithCont,holeSpace))
                P(cell,resetCell,controls(i)) = 1;
            else
                wallsCellArrival = wallsMatrix(cellArrWithCont,:);
                %apply DISTURBANCES (if impossible increase probability stay):
                %STAY
                P(cell,cellArrWithCont,controls(i)) = disturbanceSpace(3,3);
                %LEFT
                if(wallsCellArrival(3) == 0)
                    if(ismember(cellArrWithCont-M,holeSpace))
                        P(cell,resetCell,controls(i)) = disturbanceSpace(1,3); 
                    else
                        P(cell,cellArrWithCont - M,controls(i)) = disturbanceSpace(1,3); 
                    end
                else
                    P(cell,cellArrWithCont,controls(i)) = P(cell,cellArrWithCont,controls(i)) + disturbanceSpace(1,3); 
                end
                %BOTTOM
                if(wallsCellArrival(4) == 0)
                    if(ismember(cellArrWithCont-1,holeSpace))
                        P(cell,resetCell,controls(i)) = disturbanceSpace(2,3);
                    else
                        P(cell,cellArrWithCont - 1,controls(i)) = disturbanceSpace(2,3); 
                    end
                else
                    P(cell,cellArrWithCont,controls(i)) = P(cell,cellArrWithCont,controls(i)) + disturbanceSpace(2,3);  
                end
                %UP
                if(wallsCellArrival(2) == 0)
                    if(ismember(cellArrWithCont+1,holeSpace))
                        P(cell,resetCell,controls(i)) = disturbanceSpace(4,3);                        
                    else
                        P(cell,cellArrWithCont + 1,controls(i)) = disturbanceSpace(4,3);  
                    end
                else
                    P(cell,cellArrWithCont,controls(i)) = P(cell,cellArrWithCont,controls(i)) + disturbanceSpace(4,3);  
                end
                %RIGHT
                if(wallsCellArrival(1) == 0)
                    if(ismember(cellArrWithCont+M,holeSpace))
                        P(cell,resetCell,controls(i)) = disturbanceSpace(5,3); 
                    else
                        P(cell,cellArrWithCont + M,controls(i)) = disturbanceSpace(5,3);                     
                    end
                else
                    P(cell,cellArrWithCont,controls(i)) = P(cell,cellArrWithCont,controls(i)) + disturbanceSpace(5,3); 
                end
            end
        end
    end
end

end

function controls = applicableControls(cell,wallsMatrix,M,holeSpace)
%APPLICABLECONTROLS function to find feasible controls for one specific cell
%
%   Input arguments:
%
%       cell:
%           An integer containing the index of the cell that we are 
%           analyzing
%
%     wallsMatrix:
%           A (MN x 4) matrix containing, for each cell, 4 boolean values 
%          that express which walls are active. 
%          in particular:
%               1 if the wall exists, 
%               0 otherwhise
%          each value refers to a specific wall of the cell, in the 
%		   following order: 
%                        [RIGHT,UP,LEFT,BOTTOM]
%
%       M:
%           An integer containing the number of rows of the maze
%
%    	holeSpace:
%         	A (1 x H) vector containing the H holes of the maze. Each cell
%         	represents the index of a hole in the MN matrix.
%
%   Output arguments:
%
%       controls: 
%          A (C x 1) vector containing all the indexes of the avilable
%          controls for the specified cell, according to the position of
%          the walls. 
%          The indexes of the controls refer to the control space given
%          with the template
%
%	Be aware that the solution proposed here is dependent on the structure
%	of the control space provided.
%

%STAY CONTROL
controls = 7; 
%RIGHT CONTROLS
if (wallsMatrix(cell,1) == 0)
   controls = [controls, 11]; 
   cellRight = cell + M; 
   %2RIGHT
   if(wallsMatrix(cellRight,1) == 0 && not(ismember(cellRight,holeSpace)))
      controls = [controls, 13]; 
   end
end
%UP CONTROLS
if (wallsMatrix(cell,2) == 0)
   controls = [controls, 8]; 
   cellUp = cell + 1; 
   %2UP
   if(wallsMatrix(cellUp,2) == 0 && not(ismember(cellUp,holeSpace)))
      controls = [controls, 9]; 
   end
end
%LEFT CONTROLS
if (wallsMatrix(cell,3) == 0)
   controls = [controls, 3];
   cellLeft = cell - M;
   %2LEFT
   if(wallsMatrix(cellLeft,3) == 0 && not(ismember(cellLeft,holeSpace)))
      controls = [controls, 1]; 
   end
end
%BOTTOM CONTROLS
if (wallsMatrix(cell,4) == 0)
   controls = [controls, 6]; 
   cellBottom = cell - 1; 
   %2BOTTOM
   if(wallsMatrix(cellBottom,4) == 0 && not(ismember(cellBottom,holeSpace)))
      controls = [controls, 5];  
   end
end
%DIAGONALS CONTROLS RIGHT-UP
if (wallsMatrix(cell,1) == 0 && wallsMatrix(cell,2) == 0)
   cellRight = cell + M;
   cellUp = cell + 1;
   if(wallsMatrix(cellRight,2) == 0 && wallsMatrix(cellUp,1) == 0)
      controls = [controls, 12];
   end
end
%DIAGONALS CONTROLS UP-LEFT
if (wallsMatrix(cell,2) == 0 && wallsMatrix(cell,3) == 0)
   cellUp = cell + 1;
   cellLeft = cell - M;
   if(wallsMatrix(cellUp,3) == 0 && wallsMatrix(cellLeft,2) == 0)
      controls = [controls, 4];
   end
end
%DIAGONALS CONTROLS LEFT-BOTTOM
if (wallsMatrix(cell,3) == 0 && wallsMatrix(cell,4) == 0)
   cellLeft = cell - M;
   cellBottom = cell - 1;
   if(wallsMatrix(cellLeft,4) == 0 && wallsMatrix(cellBottom,3) == 0)
      controls = [controls, 2]; 
   end
end
%DIAGONALS CONTROLS BOTTOM-RIGHT
if (wallsMatrix(cell,4) == 0 && wallsMatrix(cell,1) == 0)
   cellBottom = cell - 1;
   cellRight = cell + M;
   if(wallsMatrix(cellBottom,1) == 0 && wallsMatrix(cellRight,4) == 0)
      controls = [controls, 10];
   end
end
end

function W = GenerateWallsMatrix(mazeSize, walls)
%GENERATEWALLSMATRIX function to find the active walls for each cell
%   
%	Input arguments:
%
%       mazeSize:
%           A (1 x 2) matrix containing the width and the height of the
%           maze in number of cells.
%
%       walls:
%           A (2 x 2K) matrix containing the K wall segments, where the start
%           and end point of the k-th segment are stored in column 2k-1
%           and 2k, respectively.
%
%   Output arguments:
%
%       W: 
%          A (MN x 4) matrix containing, foreach cell, 4 boolean values 
%          that express for each wall of the cell if it is active or not, 
%          in particular:
%               1 if the wall exists, 
%               0 otherwhise
%          each value refers to a specific wall of the cell, in the
%		   following order: 
%                        [RIGHT,UP,LEFT,BOTTOM]

M = mazeSize(2);
N = mazeSize(1);
K = size(walls,2)/2;
W = zeros(M*N,4);

% Set the left, right, up and bottom boundaries of the maze
for i = 1:M
    W(i,3) = 1;     
    W(M*(N-1)+i,1) = 1;
end

for i = 1:N
    W(M*(i-1)+M,2) = 1;
    W(M*(i-1)+1,4) = 1;
end

% Iterate over walls. For each wall find the cells that are contiguous to
% the wall.
for i = 1:K
    a = walls(:,2*i) - walls(:,2*i - 1);
	
    if a(1) < 0              % horizontal wall: oriented from right to left
        n = walls(1,2*i-1);
        m = walls(2,2*i-1);
        W(M*(n-1)+m,2) = 1;         % looking at cell below the wall
        W(M*(n-1)+m+1,4) = 1;       % looking at cell above the wall
		
    elseif a(1) > 0          % horizontal wall: oriented from left to right
        n = walls(1,2*i);
        m = walls(2,2*i);
        W(M*(n-1)+m,2) = 1;         % looking at cell below the wall
        W(M*(n-1)+m+1,4) = 1;       % looking at cell above the wall
		
    elseif a(2) < 0          % vertical wall: oriented  from top to bottom
        n = walls(1,2*i-1);
        m = walls(2,2*i-1);
        W(M*(n-1)+m,1) = 1;      % looking at cell on the left of the wall
        W(M*(n)+m,3) = 1;          % looking at cell on the right of the wall
		
    elseif a(2) > 0          % vertical wall: oriented  from bottom to top
        n = walls(1,2*i);
        m = walls(2,2*i);
        W(M*(n-1)+m,1) = 1;      % looking at cell on the left of the wall
        W(M*(n)+m,3) = 1;          % looking at cell on the right of the wall
    end
end
end