function outdegree_analysis(model, outdegree)
    %% Theoretical
    expected_outdegree = zeros(model.n+1,1);
    for i=2:model.n+1
        expected_outdegree(i)=expected_outdegree(i-1)+1/(i-1);
    end

    outdegree_prob_dist = zeros(model.n+1,model.n);
    for i=2:model.n+1
        outdegree_prob_dist(i,1) = 1/(i-1);
    end
    for i=3:model.n+1
        for k=2:i-1
            for l=1:i-1
                outdegree_prob_dist(i,k) = outdegree_prob_dist(i,k) + outdegree_prob_dist(l, k-1);
            end
            outdegree_prob_dist(i,k) = outdegree_prob_dist(i,k)/(i-1);
        end
    end

    exp_outdegree_squared = zeros(model.n+1,1);
    for i=1:model.n+1
        for j=1:i-1
            exp_outdegree_squared(i) = exp_outdegree_squared(i)+j.^2*outdegree_prob_dist(i,j);
        end
    end

    variance_outdegree = zeros(model.n+1,1);
    for i=1:model.n+1
        for j=1:i-1
            variance_outdegree(i) = variance_outdegree(i)+(j-expected_outdegree(i))^2*outdegree_prob_dist(i,j);
        end
    end

    outdegree_std = sqrt(variance_outdegree);

    figure();
    errorbar(expected_outdegree, outdegree_std);
    hold on;
    xlim([1,model.n-1]);
    plot(log(1:model.n+1)/log(2));
    xlabel('n');
    ylabel('Expected Outdegree');
    legend('Outdegree', 'log_2(n)');
    title('Theoretical Expected Outdegree with std');
    matlab2tikz('figures/expected_outdegree.tikz');
    
    %% Simulation
    n_sim = size(outdegree,2);
    out_resize = reshape(outdegree, [],1);
    n_links = sum(sum(outdegree(:,:)));
    average_outdegree = n_links/(model.n*n_sim);
    ratio_mean = average_outdegree/expected_outdegree(model.n);
    numerical_outdegree_variance = sqrt(1/size(out_resize,1)*ones(1,size(out_resize,1))*((out_resize-average_outdegree).^2));
    theoretical_outdegree_variance = outdegree_std(model.n-1);
    ratio_std = numerical_outdegree_variance/theoretical_outdegree_variance;
    
    figure();
    %histogram(out_resize, 0.5:model.n-0.5);
    hold on;
    %plot(1:model.n-1, outdegree_prob_dist(model.n, :)*model.n*n_sim, 'linewidth', 2);
    N = [10, 50, 100, ...%500, ...
        model.n];
    for i=1:size(N,2)
        plot(1:N(i), outdegree_prob_dist(N(i), 1:N(i)), 'linewidth', 2);
    end
    xlim([1,model.n-1]);
    xlabel('Outdegree');
    ylabel('Frequency');
    title('Outdegree distribution');
    legend('N=10', 'N=50', 'N=100', 'N=500', 'N=1000');
    matlab2tikz('figures/outdegree_pdf.tikz');
    
    
    
    


