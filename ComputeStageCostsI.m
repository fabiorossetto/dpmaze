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

% put your code here

end


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