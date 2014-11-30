function [ J_opt, u_opt_ind ] = ValueIteration( P, G )
%VALUEITERATION Value iteration
%   Solve a stochastic shortest path problem by value iteration.
%
%   [J_opt, u_opt_ind] = ValueIteration(P, G) computes the optimal cost and
%   the optimal control input for each state of the state space.
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


% Convergence tolerance (stop if ||J_k+1 - J_k|| < conv_tol)
conv_tol = 1e-4;

% Maximum number of iterations allowed
max_it = 1e4;

% Define sizes for convenience
MN = size(P,1);
L  = size(P,3);

% Initial guess for u
u = zeros(MN,1); 

J   = zeros(MN,1);
Jp1 = 10*ones(MN,1); % Silly initial guess (improve!)

it = 0;
while norm(Jp1-J) > conv_tol && it < max_it

	it = it + 1;
	J = Jp1;
	
	% This is the exact translation of the formula
	%
	% J_k+1(i) = min ( g(i,u) + sum p_ij(u) J_k(j) )
	%
	for i = 1:MN
		[Jp1(i),u(i)] = min(G(i,:) + sum(squeeze(P(i,:,:)) .* repmat(J,1,L)));
	end
	
end

J_opt = Jp1';
u_opt_ind = u';