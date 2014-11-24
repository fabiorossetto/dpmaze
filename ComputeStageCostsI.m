function G = ComputeStageCostsI( stateSpace, controlSpace, disturbanceSpace, mazeSize, walls, targetCell )
%COMPUTESTAGECOSTSI Compute stage costs.
% 	Compute the stage costs for all states in the state space for all
%   attainable control inputs.
%
%   G = ComputeStageCostsI(stateSpace, controlSpace, disturbanceSpace,
%   mazeSize, walls, targetCell) computes the stage costs for all states in
%   the state space for all attainable control inputs.
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
%       G:
%           A (MN x L) matrix containing the stage costs of all states in
%           the state space for all attainable control inputs. The entry
%           G(i, l) represents the cost if we are in state i and apply
%           control input l.

%inizialize some useful dimension
MN = size(stateSpace,1);
M = mazeSize(1);
L = size(controlSpace,1);
G = Inf(MN,L);
%create the matrix of the WALLS
wallsMatrix = GenerateWallsMatrix(mazeSize, walls);
for cell = 1:MN
    %check all AVALIAIBLE CONTROLS
    controls = availableControls(cell,wallsMatrix,M);
    for i = 1:size(controls,2)
        G(cell,controls(i)) = 1;
    end
end
%if cell is target the cost of all controls is 0
G((targetCell(2)-1)*M +targetCell(1),:) = zeros(1,L);
end

%% FUNCTION TO FIND AVAILABLE CONTROLS FOR ONE SPECIFIC CELL
%   Input arguments:
%
%       cell:
%           An integer containing the index of the cell that we are 
%           analyzing
%
%   	wallsMatrix:
%          	A (MN x 4) matrix containing, foreach cell, 4 boolean values 
%          that express foreach wall of the cell if it is active or not, 
%          in particular:
%               1 if the wall exists, 
%               0 otherwhise
%          each value referes to a specific wall of the cell, they follow 
%          this order: 
%                        [RIGHT,UP,LEFT,BOTTOM]
%
%       M:
%           An integer containing the number of rows of the maze
%
%   Output arguments:
%
%       controls: 
%          A (C x 1) vector containing all the indexes of the avilable
%          controls for the specificated cell, according to the position of
%          the walls. 
%          The indexes of the controls refere to the control space buildt
%          with the  function ScriptI
%

function controls = availableControls(cell,wallsMatrix,M)
controls = 7;
%RIGHT CONTROLS
if (wallsMatrix(cell,1) == 0)
   controls = [controls, 11];
   cellRight = cell + M; 
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
   cellLeft = cell - M; 
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

%% FUNCTION TO FIND THE ATTIVATE WALLS FOR EACH CELL
%   Input arguments:
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
%   Output arguments:
%
%       W: 
%          A (MN x 4) matrix containing, foreach cell, 4 boolean values 
%          that express foreach wall of the cell if it is active or not, 
%          in particular:
%               1 if the wall exists, 
%               0 otherwhise
%          each value referes to a specific wall of the cell, they follow 
%          this order: 
%                        [RIGHT,UP,LEFT,BOTTOM]

function W = GenerateWallsMatrix(mazeSize, walls)
M = mazeSize(1);
N = mazeSize(2);
K = size(walls,2)/2;
W = zeros(M*N,4);
for i = 1:M
    W(i,3) = 1;          % left wall
    W(M*(N-1)+i,1) = 1;   % right wall
end
for i = 1:N
    W(M*(i-1)+M,2) = 1;
    W(M*(i-1)+1,4) = 1;
end
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
    else
        disp('NEE-NOO-NEE-NOO-NEE-NOO - error: zero length wall found! NEE-NOO-NEE-NOO-NEE-NOO')
    end
end
end
