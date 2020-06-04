function [expected_clustering_coefficient, cumulative_probability] = compute_expected_clustering_coefficient(N, alpha)
expected_clustering_coefficient = zeros(N,3);
cumulative_probability = zeros(N,1);
closure_matrix = closure_probability(N);

expected_outdegree = zeros(N+1,1);
for i=2:N+1
    expected_outdegree(i)=expected_outdegree(i-1)+1/(i-1);
end

outdegree_prob_dist = zeros(N+1,N);
for i=2:N+1
    outdegree_prob_dist(i,1) = 1/(i-1);
end
cumulative_prob_dist = zeros(N+1,N);
cumulative_prob_dist(:,1) = outdegree_prob_dist(:,1);
quantile_95 = zeros(N+1,1); 
for i=3:N+1
    for k=2:i-1
        for l=1:i-1
            outdegree_prob_dist(i,k) = outdegree_prob_dist(i,k) + outdegree_prob_dist(l, k-1);
        end
        outdegree_prob_dist(i,k) = outdegree_prob_dist(i,k)/(i-1);
        cumulative_prob_dist(i,k) = cumulative_prob_dist(i,k-1) + outdegree_prob_dist(i,k);
        if (cumulative_prob_dist(i,k)>=1-alpha && quantile_95(i)==0)
            quantile_95(i)=k;
        end
    end
    for k=i:N
       cumulative_prob_dist(i,k) = cumulative_prob_dist(i,k-1) + outdegree_prob_dist(i,k);
    end
end

variance_outdegree = zeros(N+1,1);
for i=1:N+1
    for j=1:i-1
        variance_outdegree(i) = variance_outdegree(i)+(j-expected_outdegree(i))^2*outdegree_prob_dist(i,j);
    end
end

max_L = quantile_95(end);
min_p = 1/factorial(12);


tree(1).level = 1;
tree(1).p = 1;
tree(1).expected_closed_triangles = 0;
tree(1).children = [];
tree(1).parent = NaN;
tree(1).y      = 1;  
tree(1).list   = 1;
tree(1).list_size = numel(tree(1).list);
tree(1).n_triangles = 0;
tree(1).expected_clustering_coefficient = 0;

tic
[tree, ~, expected_clustering_coefficient, cumulative_probability] = add_child(tree, 1, N, 1, 0, expected_clustering_coefficient, cumulative_probability, max_L, closure_matrix, min_p);
[tree, ~, expected_clustering_coefficient, cumulative_probability] = add_child(tree, 1, N, 1, 1, expected_clustering_coefficient, cumulative_probability, max_L, closure_matrix, min_p);
toc

% theoretical_triangles = (expected_outdegree.*(expected_outdegree-1) + variance_outdegree)/2;
% 
% figure();
% plot(2:N+1, expected_clustering_coefficient(:,1))
% hold on;
% plot(2:N+1, expected_clustering_coefficient(:,2))
% plot(2:N+1, theoretical_triangles(2:N+1), '+');
% legend('expected triangles', 'expected closed triangles', 'theoretical expected triangles');
% legend('expected triangles', 'expected closed triangles', 'theoretical expected triangles');
% hold off;

