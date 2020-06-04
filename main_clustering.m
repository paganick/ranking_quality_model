clean all;
close all;
clc;

N = 30;
n_sim = 10;
average_clustering = NaN(N+1,1);
alpha = 0.05;
nn  = 3:5:N+1;
nn  = 1000;
%[expected_clustering_coefficient, cumulative_probability] = compute_expected_clustering_coefficient(N, alpha);

for i=1:size(nn,2)
    model.n = nn(i);
    model.T = 2000;
    model.q      = rand(model.n,1);
    model.q      = sort(model.q, 'descend');

    sim_average_clustering =  zeros(n_sim, N+1);
    tic
    for sim=1:n_sim
        A = simulate(model);
        clustering_coefficient = compute_clustering(A);
        sim_average_clustering(sim, nn(i)) = mean(clustering_coefficient);
        %sim_average_clustering(sim, nn) = (clustering_coefficient(end));
        scatter(1:model.n, clustering_coefficient);
    end
    toc
    average_clustering(nn(i)) = mean(sim_average_clustering(:, nn(i)));
    %scatter(3:N+1, average_clustering(3:end));
    %pause(0.01);
end







%%
% errorpos = (1-cumulative_probability);
% errorneg = zeros(N+1,1);
% figure();
% errorbar(3:N+1, expected_clustering_coefficient(2:end,3), errorneg(3:N+1), errorpos(2:end), 'red');
% hold on;
% scatter(3:N+1, average_clustering(3:end));
% legend('Expected theoretical', 'Numerical average');
% title('Average Clustering');
% nn = [3:1:N+1]';
% sup_estimate = expected_clustering_coefficient(2:end,3)+errorpos(2:end);
% output_data = [nn, expected_clustering_coefficient(2:end,3), sup_estimate];
% save('result.txt', 'output_data', '-ascii')

%%
mean_value = expected_clustering_coefficient(2:end,3)./(cumulative_probability(2:end));
errorpos = (1-cumulative_probability)+ expected_clustering_coefficient(:,3).*(1-1./(cumulative_probability));
errorneg = expected_clustering_coefficient(2:end,3) - mean_value;
figure();
errorbar(3:N+1, mean_value, errorneg, errorpos(2:end), 'red');
hold on;
scatter(3:N+1, average_clustering(3:end));
legend('Expected theoretical', 'Numerical average');
title('Average Clustering');
nn = [3:1:N+1]';
output_data2 = [nn, mean_value, errorpos(2:end), errorneg];
save('result2.txt', 'output_data2', '-ascii');

%%
% derivative = zeros(N-1, 1);
% sec_derivative = zeros(N-2, 1);
% for i=1:N-1
%     derivative(i) = expected_clustering_coefficient(i+1,3)-expected_clustering_coefficient(i,3);
% end
% for i=1:N-2
%     sec_derivative(i) = derivative(i+1)-derivative(i);
% end
% figure();
% plot(2:N, derivative);
% title('First derivative');
% figure();
% plot(3:N, sec_derivative);
% title('Second derivative');
