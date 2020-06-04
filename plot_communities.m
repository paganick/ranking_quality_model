close all;
clear all;
clc;
T_comm = readtable("data/nodes_with_communities.csv");
T = readtable("data/edges.csv");
s = table2array(T(:,1));
t = table2array(T(:,2));
community = table2array(T_comm(:,4));
s = s+1;
t = t+1;
community = community+1;
n_edges = size(s,1);
communities = unique(community);
n_communities = size(communities,1);
[size_communities, ~] = hist(community, n_communities-1);
ss = cell(1, n_communities);
tt = cell(1, n_communities);
for i=1:n_edges
    if ( community(s(i)) == community(t(i)) ) %extreme nodes are in the same community
        ss{community(s(i))} = [ss{community(s(i))}; s(i)];
        tt{community(t(i))} = [tt{community(t(i))}; t(i)];
    end
end

for i=1:n_communities
    A = [];
    figure(2*(i-1)+1);
    nodes = unique([ss{i};tt{i}]);
    n = size(nodes,1)
    G=digraph( ss{i}, tt{i} );
    A = adjacency(G);
    %plot(G);
    %A( ~any(A,2), : ) = [];  %rows
    %A( :, ~any(A,1) ) = [];  %columns
    %A = full(A);
    %spy(A);
    indegree = sum(A,1);
    [max_indeg, max_index] = max(full(indegree));
    hub_connection_index = find(A(:, max_index));
    A(~hub_connection_index, :) = [];
    A(:, ~hub_connection_index) = [];
    n = size(hub_connection_index,1);
    outdegree = sum(A,2);
    exp_outdegree = harmonic(n)
    average_outdegree = size(ss{i},1)/n
    indegree = sum(A,1);
    indegree = sort(indegree, 'descend');
    %figure(2*i);
    scatter(1:min(20,n), indegree(1:min(20,n)));
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
end


