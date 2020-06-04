clear all;
clean all;
close all;
%clc;


n = 100;

connection_probability = zeros(n,n);

for i=1:n
    for j=1:i-1
        connection_probability(i,j) = 1/i + 1/j - 1/(i*j);
    end
end

connection_probability = connection_probability+connection_probability';


expected_clustering_coefficient = zeros(n,3);
cumulative_probability = zeros(n,1);

expected_outdegree = zeros(n+1,1);
for i=2:n+1
    expected_outdegree(i)=expected_outdegree(i-1)+1/(i-1);
end

outdegree_prob_dist = zeros(n+1,n);
for i=2:n+1
    outdegree_prob_dist(i,1) = 1/(i-1);
end
cumulative_prob_dist = zeros(n+1,n);
cumulative_prob_dist(:,1) = outdegree_prob_dist(:,1);
alpha = 0.05;
quantile_95 = zeros(n+1,1); 
for i=3:n+1
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
    for k=i:n
       cumulative_prob_dist(i,k) = cumulative_prob_dist(i,k-1) + outdegree_prob_dist(i,k);
    end
end

max_L = quantile_95(end);
min_p = 1/factorial(12);

variance_outdegree = zeros(n+1,1);
for i=1:n+1
    for j=1:i-1
        variance_outdegree(i) = variance_outdegree(i)+(j-expected_outdegree(i))^2*outdegree_prob_dist(i,j);
    end
end

outdegree_std = sqrt(variance_outdegree);

figure();
errorbar(expected_outdegree, outdegree_std);
hold on;
plot(log([1:n+1])/log(2));
legend('Outdegree', 'log_2(n)');
title('Expected Outdegree');

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
[tree, ~, expected_clustering_coefficient, cumulative_probability] = add_child(tree, 1, n, 1, 0, expected_clustering_coefficient, cumulative_probability, max_L, connection_probability, min_p);
[tree, ~, expected_clustering_coefficient, cumulative_probability] = add_child(tree, 1, n, 1, 1, expected_clustering_coefficient, cumulative_probability, max_L, connection_probability, min_p);
toc


theoretical_triangles = (expected_outdegree.*(expected_outdegree-1) + variance_outdegree)/2;

% figure();
% plot(2:n+1, expected_clustering_coefficient(:,1))
% hold on;
% plot(2:n+1, expected_clustering_coefficient(:,2))
% plot(2:n+1, theoretical_triangles(2:n+1), '+');
% legend('expected triangles', 'expected closed triangles', 'theoretical expected triangles');
% legend('expected triangles', 'expected closed triangles', 'theoretical expected triangles');
% hold off;

%errorpos = (1-cumulative_prob_dist(:,max_L));
errorpos = (1-cumulative_probability);
errorneg = zeros(n+1,1);
figure();
errorbar(3:n+1, expected_clustering_coefficient(2:end,3), errorneg(3:n+1), errorpos(2:end), 'red');
hold on;
%approx_clustering_coefficient = expected_clustering_coefficient(:,2)./expected_clustering_coefficient(:,1);
%plot(3:n+1, approx_clustering_coefficient(2:end), 'b');
hold off;
title('Clustering');
%legend('theoretical', 'approx');
%legend('theoretical', 'approx', '1^{st} derivative', '2^{nd} derivative');

derivative = zeros(n-1, 1);
sec_derivative = zeros(n-2, 1);
for i=1:n-1
    derivative(i) = expected_clustering_coefficient(i+1,3)-expected_clustering_coefficient(i,3);
end
for i=1:n-2
    sec_derivative(i) = derivative(i+1)-derivative(i);
end
% plot(2:n, derivative, 'g+');
% plot(3:n, sec_derivative, 'magenta');


figure();
plot(2:n, derivative);

figure();
plot(3:n, sec_derivative);

% x = [2:n]';
% y = approx_clustering_coefficient(2:end);
% x0 = [1 1 1];
% %fitfun = fittype( @(a,b,c, x) a*exp(b*(x))+c);
% fitfun = fittype( @(a,b,c, x) a./(x+b)+c);
% [fitted_curve,gof] = fit(x,y,fitfun,'StartPoint',x0)
% figure()
% xx = [0:100]';
% plot(x, y, 'p',  xx, fitted_curve(xx), '-r')
% %plot(x, y, 'p',  x, a*exp(b*(x-2)), '-r')

nn = [3:1:n+1]';
sup_estimate = expected_clustering_coefficient(2:end,3)+errorpos(2:end);
output_data = [nn, expected_clustering_coefficient(2:end,3), sup_estimate];
save('result.txt', 'output_data', '-ascii')


mean_value = expected_clustering_coefficient(2:end,3)./(cumulative_probability(2:end));
errorpos = (1-cumulative_probability)+ expected_clustering_coefficient(:,3).*(1-1./(cumulative_probability));
errorneg = expected_clustering_coefficient(2:end,3) - mean_value;
figure();
errorbar(3:n+1, mean_value, errorneg, errorpos(2:end), 'red');
hold on;
%approx_clustering_coefficient = expected_clustering_coefficient(:,2)./expected_clustering_coefficient(:,1);
%plot(3:n+1, approx_clustering_coefficient(2:end), 'b');
hold off;
title('Clustering');

output_data2 = [nn, mean_value, errorpos(2:end), errorneg];
save('result2.txt', 'output_data2', '-ascii');
