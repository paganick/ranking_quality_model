clear all;
close all;
clc;

model.n = 1000;

% The outdegree is set up similar to the average outdegree of the Ranking
% quality model.
model.outdegree = ceil(harmonic(model.n));
n_sim = 10;

final_indegree = zeros(model.n, n_sim);
for sim=1:n_sim
    A = SF_simulate(model);
    indegree_sorted = sort(sum(A,1), 'descend')';
    final_indegree(:, sim) = indegree_sorted;
end
final_indegree = final_indegree - model.outdegree; % Removing the initial
% bi-lateral connections which are an artifact.
SF_indegree_analysis(model, final_indegree);

