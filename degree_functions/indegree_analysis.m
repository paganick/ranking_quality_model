function indegree_analysis(model, indegree)

    n_sim = size(indegree,2);
    indegree_prob_table = compute_indegree_prob_table(model);
    theoretical_indegree_pdf = sum(indegree_prob_table,1)/model.n;
    theoretical_indegree_cdf = cumsum(theoretical_indegree_pdf);
    
    numerical_pdf = zeros(model.n,1);
    for d=1:model.n
        numerical_pdf(d) = histcounts(indegree(:,:), [d-1 d-1])/(model.n*n_sim);
    end
    numerical_indegree_cdf  = cumsum(numerical_pdf);

    %%
    % Plot 1: Indegree Probability Distribution Function, Theoretical and
    % Numerical.
    figure();
    hold on;
    x = 0:1:model.n-1;
    scatter(x, theoretical_indegree_pdf,'+');
    scatter(x, numerical_pdf, 'd');
    power_law_pdf_fit = fit(log(1:(model.n-1)/2+1)', log(theoretical_indegree_pdf(1:(model.n-1)/2+1)'), 'poly1');
    y=[];
    for d=2:model.n-1
        y = [y; repmat([d-1], ceil(1000*theoretical_indegree_pdf(d+1)),1)];
    end
    lognfit_pdf_fit = fitdist(y+1, 'lognormal'); %TO BE CHECKED
    xlabel('Indegree');
    ylabel('pdf');
    xlim([0,model.n-1]);
    plot(x, lognpdf(x, lognfit_pdf_fit.mu-1, lognfit_pdf_fit.sigma)); %TO BE CHECKED
    plot(1:1:model.n-1, exp(power_law_pdf_fit(log(1:1:model.n-1))));
    log_normal_fit_legend = ['log-normal fit: \mu=', num2str(lognfit_pdf_fit.mu-1), ', \sigma=', num2str(lognfit_pdf_fit.sigma)];
    power_law_fit_legend = ['power-law fit: p1=', num2str(power_law_pdf_fit.p1), ', p2=', num2str(power_law_pdf_fit.p2)];
    % p1 should be close to -2, for pure Zipf's law, or for a power law
    % with coeff. -2.
    legend('Theoretical at T=\infty', 'Numerical', log_normal_fit_legend, power_law_fit_legend);
    title('Indegree Pdf');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');

    %% 
    % Plot 2: Indegree Survival Distribution Function, Theoretical and
    % Numerical.
    
    figure();
    hold on;
    stairs(0:model.n-1, 1-numerical_indegree_cdf, 'LineWidth', 2);
    stairs(0:model.n-1, 1-theoretical_indegree_cdf, 'LineWidth', 2);
    set(gca, 'XScale', 'log');
    xlabel('Indegree');
    ylabel('Survival Dist. Funct.');
    legend('Empirical', 'Theoretical at T=\infty');
    xlim([0,model.n-1]);
    set(gca,'Ylim',[10^-2, 1]);
    title('Indegree Survival Distribution Function');
    set(gca, 'YScale', 'log');

    %% 
    % Plot: Indegree Probability Distribution Function vs Rank
    % According to Zipf's law, expected slope of -1.
    figure();
    [D,K]=meshgrid(1:1:model.n, 0:1:model.n-1);
    surf(D,K,indegree_prob_table');
    shading interp;
    colormap(flipud(hot));
    hcb=colorbar();
    hcb.Label.String = 'Probability';
    hold on;
    
    average_indegree = mean(indegree, 2);
    scatter3(1:model.n, average_indegree, ones(1,model.n), 'b', 'filled');

    for i=1:model.n
        expected_mean_indegree(i) = exp_mean_indegree(model, i);
    end
    plot3(1:model.n, expected_mean_indegree, ones(1,model.n), 'g', 'LineWidth', 2);
    
    expected_mean_indegree_inf= (model.n+1)*ones(model.n,1) - [model.n:-1:1]';
    expected_mean_indegree_inf= model.n./expected_mean_indegree_inf;
    expected_mean_indegree_inf(1) = model.n;
    plot3(expected_mean_indegree_inf, 1:model.n, ones(1,model.n), 'm+', 'LineWidth',1);
    
    view(2);
    xlim([1,model.n]);
    ylim([1,model.n-1]);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('Rank');
    ylabel('Indegree');
    legend('Pdf', 'Numerical Mean', 'Expected Mean', 'Expected Mean at T=\infty');
    title('Zipf`s law'); 
    hold off;

    %% 
    % Plot: Fit of Zipf's law from Numerical Data.
    n_min = 1;
    n_max = find(~average_indegree, 1, 'first')-1;
    if (isempty(n_max))
        n_max = model.n;
    end
    [f1, f2] = fit(log(n_min:n_max)', log(average_indegree(n_min:n_max)), 'poly1');

    figure();
    scatter(1:model.n, average_indegree(1:model.n), 'b', 'filled');
    hold on;
    plot(1:model.n, expected_mean_indegree(1:model.n), 'g', 'LineWidth',2);
    plot(n_min:n_max, exp(f1(log(n_min:n_max))), 'LineWidth', 2);
    xlabel('Rank');
    ylabel('Indegree');
    xlim([0,model.n-1]);
    ylim([0,model.n-1]);
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    power_law_fit_legend = ['Power-law fit: p1=', num2str(f1.p1), ', p2=', num2str(f1.p2)];
    legend('Numerical Mean', 'Expected Mean', power_law_fit_legend);
    title('Power-law fit of Zipf`s law');
    
