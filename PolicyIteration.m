function [ J_opt, u_opt_ind ] = PolicyIteration( P, G )
%POLICYITERATION Value iteration
%   Solve a stochastic shortest path problem by policy iteration.
%
%   [J_opt, u_opt_ind] = PolicyIteration(P, G) computes the optimal cost and
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

% Initial guess for mu. Greedy strategy: choose the locally cheaper control
mu = 7*ones(MN-1,1); % min(G,[],2);
mu_m1 = 10*ones(MN-1,1); 

[t,~] = find(G == 0,1); % Terminal state

Pmt = P([1:t-1 , t+1:end],[1:t-1 , t+1:end],:);
Gmt = G([1:t-1 , t+1:end],:);

it = 0;
while norm(mu_m1-mu) > conv_tol && it < max_it

	it = it + 1
	mu_m1 = mu;
	
	% Solve linear system
	J = (eye(MN-1)-prob(Pmt,mu)) \ cost(Gmt,mu);
	
	% Improve
	for i = 1:MN-1
		[~,mu(i)] = min(Gmt(i,:) + sum(squeeze(Pmt(i,:,:)) .* repmat(J,1,L)));
	end
end

J_opt = [J(1:t-1) ; 0 ; J(t:end)]';
u_opt_ind = [mu(1:t-1) ; 7 ; mu(t:end)]'; %TODO remove 7
end

function Pmu = prob(P,mu)
MNm1 = length(mu);

Pmu = zeros(MNm1);
for i = 1:MNm1;
	for j = 1:MNm1;
		Pmu(i,j) = P(i,j,mu(i));
	end
end
end

function Gmu = cost(G,mu)
MNm1 = length(mu);
Gmu = zeros(MNm1,1);
for i = 1:MNm1
	Gmu(i) = G(i,mu(i));
end
end
