function SF_indegree_analysis(model, final_indegree)
    
    %% Probability Distribtion Function plot.
    n_sim = size(final_indegree,2);
    for d=1:model.n
        numerical_pdf(d) = histcounts(final_indegree(:,:), [d-1 d-1])/(model.n*n_sim);
    end
    figure();
    hold on;
    x = 0:1:model.n-1;
    scatter(x, numerical_pdf, 'filled');
    [~, yy, v] = find(numerical_pdf);
    power_law_pdf_fit = fit(log(yy(2:end)-1)', log(v(2:end))', 'poly1');
    plot(yy, exp(power_law_pdf_fit(log(yy))), 'LineWidth', 2);
    xlabel('Indegree');
    ylabel('pdf');
    xlim([0,model.n-1]);
    power_law_fit_legend = ['power-law fit: p1=', num2str(power_law_pdf_fit.p1), ', p2=', num2str(power_law_pdf_fit.p2)];
    legend('Numerical', power_law_fit_legend);
    title('Indegree Pdf');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');


    %% Comparison with Zipf's law
    average_indegree = sum(final_indegree,2)/n_sim;
    % Zipf's law
    zipf_indegree= (model.n+1)*ones(model.n,1) - [model.n:-1:1]';
    zipf_indegree= model.n./zipf_indegree;
    zipf_indegree(1) = model.n;


    correlation = corr(average_indegree, zipf_indegree);

    n_min = 1;
    n_max = model.n;
    if (any(~average_indegree))
        n_max = find(~average_indegree, 1, 'first')-1;
    end

    [f1, ~] = fit(log(n_min:n_max)', log(average_indegree(n_min:n_max)), 'poly1');

    figure();
    scatter(1:n_max, average_indegree(1:n_max), 'filled');
    hold on;
    plot(n_min:n_max, exp(f1(log(n_min:n_max))), 'LineWidth', 2);
    plot(1:n_max, zipf_indegree(1:n_max), 'LineWidth', 2);
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlabel('Node rank');
    ylabel('Indegree');
    power_law_fit_legend = ['Power law fit: p1=', num2str(f1.p1), ' p2=', num2str(f1.p2)];
    legend('Numerical average', power_law_fit_legend, 'Zipf`s law (p1=-1)');
    title('Directed preferential attachment and Zipf`s law');
    hold off;

