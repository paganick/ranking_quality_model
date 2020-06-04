clear all;
close all;
clc;

N = 7;

connection_probability = zeros(N,N);

for n_i=1:N
    for j=1:n_i-1
        connection_probability(n_i,j) = 1/n_i + 1/j - 1/(n_i*j);
    end
end

connection_probability = connection_probability+connection_probability';

max_n = factorial(N);
combinations={};
probability = [];
cum_probability = [];
n_combinations = 0;
cum_n_combinations = [1;2];
cum_probability =    [1/N;2/N];
tic
c = expected_clustering([1, 2], connection_probability)*1/N;

for n=3:max_n+1
    combinations{n} = {};
    n_i=0;
    probability(n) = 1/N*(1/(n-1));
    if (n<=N)
        n_i = n_i+1;
        combinations{n}{n_i} = [n];
        c = c+expected_clustering([1, n], connection_probability)*probability(n);
        c = c+expected_clustering([1, 2, n], connection_probability)*probability(n);
    end
    div = divisors(n-1);
    div = div(div>1 & div<n-1);
    for i=1:size(div,2)
        j = div(i);
        k = (n-1)/j;
        if (k>1 && k~=j)
            sub_combinations = combinations{j+1};
            for i_comb=1:size(sub_combinations,2)
                new_combination = [sub_combinations{i_comb}, k+1];
                if (combination_is_good(new_combination, N))
                    n_i = n_i+1;
                    combinations{n}{n_i} = new_combination;
                    c = c+expected_clustering([1, new_combination], connection_probability)*probability(n);
                    c = c+expected_clustering([1, 2, new_combination], connection_probability)*probability(n);
                end
            end
        end
    end
    n_combinations(n) = size(combinations{n},2)*2;
    cum_n_combinations(n) = cum_n_combinations(n-1) + n_combinations(n);
    cum_probability(n)  = cum_probability(n-1)+probability(n)*n_combinations(n);
    if (cum_probability(n)>=1.00)
        break;
    end
end
toc

c

figure();
plot(cum_probability);

function c = expected_clustering(combination, connection_probability)
    L = size(combination,2);
    c = 0;
    if (L>1)
        expected_triangles = L*(L-1)/2;
        expected_closed_triangles = 0;
        for i=1:L
            for j=1:i-1
                expected_closed_triangles = expected_closed_triangles + connection_probability(combination(i),combination(j));
            end
        end
        c = expected_closed_triangles/expected_triangles;
    end
end

function good = combination_is_good(combination, N)
    if (combination(1) <= N)
        good = true;
    else
        good = false;
    end
    if (size(combination,2)>1)
        for i=2:size(combination,2)
            if (combination(i)<=combination(i-1) || combination(i)>N)
                good = false;
            end
        end
    end
end
%%

% combinations{1} = 1;
% combinations{2} = [1,2];
% L=1;
% max_L=1;
% %ii = NaN(N,1);
% ii = 3;
% p_length = zeros(N-2,1);
% for j=1:N-2
%     p_length(j) = factorial(j+2)/2;
% end
% 
% next_combination = {};
% for i=1:N-2
%     next_combination{i} = [];
% end
% for i=1:N-2
%     for j=1:i
%         next_combination{i} = [next_combination{i}, 2+j];
%     end
% end
% 
% i=3;
% while (L<=N-2)
%     [m, L] = min(p_length);
%     combinations{i} = next_combination{L}; i = i+1;
%     [next_combination{L}, found] = get_next_combination(L, next_combination{L}, N);
%     p_length(L) = compute_p(next_combination{L});
%     if (found == false)
%         p_length(L) = Inf;
%     end
%     combinations{i-1}
% end
% 
% function p = compute_p(combination)
%     p=1;
%     for i=1:size(combination,2)
%         p = p*combination(i);
%     end
% end
%     
% 
% function [next_combination, found] = get_next_combination(L, combination, N)
%     found = false;
%     next_combination = combination;
%     for i=1:L
%         if (combination(end-i+1)<N)
%             combination
%             next_combination(end-i+1) = next_combination(end-i+1)+1;
%             next_combination
%             found = true;
%             break;
%         end
%     end
% end

%% 
% exp_ci = zeros(N,1);
% p_ci   = zeros(N,1);
% 
% alpha = 0.05;
% 
% combinations={};
% combinations{1} = 1;
% combinations{2} = [1,2];
% combinations_probability = [1/2, 1/2];
% combinations_probability_factor = [1, 1];
% exp_closed_triangles(1) = connection_probability(1,1);
% exp_closed_triangles(2) = connection_probability(1,2);
% combination_length = [1, 2];
% exp_triangles = 1/2*(combination_length).*(combination_length-1);
% 
% other_combinations = {};
% 
% for n=3:N
%     new_combinations = combinations;
%     new_combinations_probability = combinations_probability*1/n;
%     new_combination_length = combination_length+1;
%     new_exp_triangles = 1/2*(new_combination_length).*(new_combination_length-1);
%     new_combinations_probability_factor = combinations_probability_factor*n;
%     new_exp_closed_triangles = exp_closed_triangles;
%     for i=1:size(combinations,2)
%         for j=1:size(new_combinations{i})
%             new_exp_closed_triangles(i) = new_exp_closed_triangles(i)+connection_probability(new_combinations{i}(j),n);
%         end
%         new_combinations{i} = [new_combinations{i}, n];
%     end
%     combinations_probability = combinations_probability*(n-1)/n;
%     new_combos_index = 1;
%     old_combos_index = 1;
%     for i=1:size(combinations,2)*2
%         old_combos_weight = combinations_probability_factor(old_combos_index);
%         new_combos_weight = new_combinations_probability_factor(new_combos_index);
%         if (old_combos_weight>new_combos_weight || new_combos_index==size(combinations,2)) 
%             p_ci(n) = p_ci(n)+ combinations_probability(new_combos_index);
%             exp_ci(n) = exp_ci(n) + exp_closed_triangles(new_combos_index)./exp_triangles(new_combos_index);
%             if (p_ci(n) > 1-alpha)
%                 break;
%             end
%             new_combos_index = new_combos_index+1;
%         else
%             p_ci(n) = p_ci(n)+ combinations_probability(old_combos_index);
%             exp_ci(n) = exp_ci(n) + exp_closed_triangles(old_combos_index)./exp_triangles(old_combos_index);
%             if (p_ci(n) > 1-alpha)
%                 break;
%             end
%             old_combos_index = old_combos_index+1;         
%         end
%     end
% end
