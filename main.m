clear all;
close all;
clc;


model.n = 201;

model.T = 10000;
n_sim = 1000;

% Quality definition
model.q      = rand(model.n,1);
model.q      = sort(model.q, 'descend');

final_indegree=zeros(model.n, n_sim);
final_outdegree=zeros(model.n, n_sim);
overlap_matrix = zeros(model.n, model.n);

tic
for sim=1:n_sim
    A = simulate(model);
    final_indegree(:, sim) = sum(A,1)';
    final_outdegree(:, sim) = sum(A,2)';
    %overlap_matrix = overlap_analysis(model, A, overlap_matrix, sim);
end
toc

%% Indegree Analysis

indegree_analysis(model, final_indegree);



%% Outdegree Analysis

outdegree_analysis(model, final_outdegree);



