function [ J_opt, u_opt_ind ] = LinearProgramming( P, G )
%LINEARPROGRAMMING Value iteration
%   Solve a stochastic shortest path problem by linear programming.
%
%   [J_opt, u_opt_ind] = LinearProgramming(P, G) computes the optimal cost
%   and the optimal control input for each state of the state space.
%
%   Input arguments:
%
%       P:
%           A (MN x MN x L) matrix containing the transition probabilities
%           between all states in the state space for all attainable
%           control inputs. The entry P(i, j, l) represents the transition
%           probability from state i to state j if control input l is
%           applied.
%
%       G:
%           A (MN x L) matrix containing the stage costs of all states in
%           the state space for all attainable control inputs. The entry
%           G(i, l) represents the cost if we are in state i and apply
%           control input l.
%
%   Output arguments:
%
%       J_opt:
%       	A (1 x MN) matrix containing the optimal cost-to-go for each
%       	element of the state space.
%
%       u_opt_ind:
%       	A (1 x MN) matrix containing the indices of the optimal control
%       	inputs for each element of the state space.

% Define sizes for convenience
MN = size(P,1);
L  = size(P,3);

[t,~] = find(G == 0,1); % Terminal state

Pmt = P([1:t-1 , t+1:end],[1:t-1 , t+1:end],:);
Gmt = G([1:t-1 , t+1:end],:);

f = -1 * ones(MN-1,1);
A = repmat(eye(MN-1),L,1) - reshape(permute(Pmt,[1,3,2]),(MN-1)*L,MN-1);
b = Gmt(:);

J = linprog(f,A,b);

u = zeros(MN-1);
for i = 1:MN-1
	[~,u(i)] = min(Gmt(i,:) + sum(squeeze(Pmt(i,:,:)) .* repmat(J,1,L)));
end

J_opt = [J(1:t-1) ; 0 ; J(t:end)]';
u_opt_ind = [u(1:t-1) ; 7 ; u(t:end)]'; %TODO remove 7

end

