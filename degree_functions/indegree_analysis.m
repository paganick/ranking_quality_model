function indegree_analysis(model, indegree)
    
    %% Plot 1:
    % fitting the CDF with the method specified in Clauset 2007.
    indegree_1d_withzero = reshape(indegree, 1, []);
    indegree_1d = indegree_1d_withzero(indegree_1d_withzero~=0);
    [alpha, xmin, L] = plfit(indegree_1d)
    [h, powerlawfit_cdf] = plplot(indegree_1d, xmin, alpha);
    [alpha_tol, xmin_tol, n_tol] = plvar(indegree_1d, 'silent')
    [p,gof]=plpva(indegree_1d, xmin, 'silent')
    power_law_fit_legend = ['power-law fit: alpha=', num2str(alpha),  ', alpha_{tol}=', num2str(alpha_tol), ', xmin=', num2str(xmin) ', xmin_{tol}=', num2str(xmin_tol),  ', p-value=', num2str(p), ', gof=', num2str(gof)];
    legend('Numerical', power_law_fit_legend);
    %matlab2tikz('figures/indegree_cdf_new.tikz');
 
    %% Fit the numerical data with lognormal distribution
    lognormal_parameter_hat = lognfit(indegree_1d);
    
    %% Computing numerical pdf and cdf
    n_sim = size(indegree,2);
    numerical_pdf = zeros(model.n,1);
    for d=1:model.n
        numerical_pdf(d) = histcounts(indegree(:,:), [d-1 d-1])/(model.n*n_sim);
    end
    numerical_cdf = cumsum(numerical_pdf);
    
    %% Computing theoretical pdf and cdf
    [indegree_prob_table, indegree_prob_table_inf] = compute_indegree_prob_table(model);
    theoretical_indegree_pdf = sum(indegree_prob_table_inf,1)/model.n;
    theoretical_indegree_cdf = cumsum(theoretical_indegree_pdf);
    
    
    %% Plot 2:
    % Numerical pdf and lognormal fit.
    figure();
    hold on;
    x = 0:model.n-1;
    %theoretical
    scatter(x, theoretical_indegree_pdf,'+');
    %numerical
    scatter(x, numerical_pdf, 'd');
    %lognormal
    pd = makedist('Lognormal','mu',lognormal_parameter_hat(1),'sigma',lognormal_parameter_hat(2));
    y = pdf(pd, x);
    plot(x,y);
    %powerlaw fit
    powerlawfit_pdf = powerlawfit_cdf;
    for i=1:size(powerlawfit_pdf,1)-1
        powerlawfit_pdf(i,2) = powerlawfit_cdf(i,2) -  powerlawfit_cdf(i+1,2);
    end
    plot(powerlawfit_pdf(:,1), powerlawfit_pdf(:,2));
    
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlim([0,model.n-1]);
    ylim([10^-6, 1]);
  
    log_normal_fit_legend = ['log-normal fit: \mu=', num2str(lognormal_parameter_hat(1)), ', \sigma=', num2str(lognormal_parameter_hat(2))];
    power_law_fit_legend = ['power-law fit: alpha=', num2str(alpha),  ', alpha_{tol}=', num2str(alpha_tol), ', xmin=', num2str(xmin) ', xmin_{tol}=', num2str(xmin_tol),  ', p-value=', num2str(p), ', gof=', num2str(gof)];
    legend('Theoretical at T=\infty', 'Numerical distribution', ...
            log_normal_fit_legend, power_law_fit_legend);
    matlab2tikz('figures/indegree_pdf.tikz');
    
    %%
    % Plot 3: Indegree Complementary Cumulative Distribution Function, 
    % Theoretical, Numerical, and fits.
    figure();
    hold on;
    x = 1:1:model.n-1;
    %theoretical
    scatter(x, 1-theoretical_indegree_cdf(1:end-1),'+');
    %numerical
    scatter(x, 1-numerical_cdf(1:end-1), 'd');
    %lognormal
    pd = makedist('Lognormal','mu',lognormal_parameter_hat(1),'sigma',lognormal_parameter_hat(2));
    y = cdf(pd, x);
    plot(x,1-y);
    %powerlaw
    plot(powerlawfit_cdf(:,1), powerlawfit_cdf(:,2));
    
    xlabel('Indegree');
    ylabel('cdf');
    xlim([0,model.n-1]);
    ylim([10^-3, 1]);
    log_normal_fit_legend = ['log-normal fit: \mu=', num2str(lognormal_parameter_hat(1)), ', \sigma=', num2str(lognormal_parameter_hat(2))];
    power_law_fit_legend = ['power-law fit: alpha=', num2str(alpha),  ', alpha_{tol}=', num2str(alpha_tol), ', xmin=', num2str(xmin) ', xmin_{tol}=', num2str(xmin_tol),  ', p-value=', num2str(p), ', gof=', num2str(gof)];
    legend('Theoretical at T=\infty', 'Numerical distribution', ...
            log_normal_fit_legend, power_law_fit_legend);
    title('Indegree CCdf');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    matlab2tikz('figures/indegree_cdf.tikz');


    %% 
    % Plot 4: Indegree Probability Distribution Function vs Rank, at time t.
    figure();
    [D,K]=meshgrid(1:1:model.n, 1:1:model.n-1);
    h = surf(D,K,indegree_prob_table(:,2:end)');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    shading interp;
    cMap = flipud(gray);
    dataMax = 1;
    dataMin = 0;
    centerPoint = 0;
    scalingIntensity = 5;
    x = 1:length(cMap); 
    x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
    x = scalingIntensity * x/max(abs(x));
    x = sign(x).* exp(abs(x));
    x = x - min(x); x = x*511/max(x)+1; 
    newMap = interp1(x, cMap, 1:512);
    colormap(newMap)
    view(2);
    xlim([1,model.n]);
    ylim([1,model.n-1]);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    set(gca,'visible','off');
    box off;
    H = getframe(gca);
    imwrite(H.cdata, 'figures/indegree_mesh_T=200.png')
    
    figure();
    expected_mean_indegree = indegree_prob_table*[0:model.n-1]';
    plot3(1:model.n, expected_mean_indegree,  ones(1,model.n), '+', 'LineWidth', 1);
    view(2);
    xlim([1,model.n]);
    ylim([1,model.n-1]);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('Rank i');
    ylabel('Indegree d_{in}^i');
    matlab2tikz('figures/mean_indegree_zipf_T=200.tikz');
    
    %% 
    % Plot 5: Indegree Probability Distribution Function vs Rank
    % According to Zipf's law, expected slope of -1.
    figure();
    [D,K]=meshgrid(1:1:model.n, 1:1:model.n-1);
    h = surf(D,K,indegree_prob_table_inf(:,2:end)');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    shading interp;
    cMap = flipud(gray);
    dataMax = 1;
    dataMin = 0;
    centerPoint = 0;
    scalingIntensity = 5;
    x = 1:length(cMap); 
    x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
    x = scalingIntensity * x/max(abs(x));
    x = sign(x).* exp(abs(x));
    x = x - min(x); x = x*511/max(x)+1; 
    newMap = interp1(x, cMap, 1:512);
    colormap(newMap);
    view(2);
    xlim([1,model.n]);
    ylim([1,model.n-1]);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    colorbar
    caxis([0 1])
    set(gca,'visible','off');
    box off;
    H = getframe(gca);
    imwrite(H.cdata, 'figures/indegree_mesh_T=inf.png')
    
    figure();
    expected_mean_indegree_inf= (model.n+1)*ones(model.n,1) - [model.n:-1:1]';
    expected_mean_indegree_inf= model.n./expected_mean_indegree_inf;
    expected_mean_indegree_inf(1) = model.n-1;
    plot3(1:model.n, expected_mean_indegree_inf,  ones(1,model.n), '+', 'LineWidth', 1);
    view(2);
    xlim([1,model.n]);
    ylim([1,model.n-1]);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('Rank i');
    ylabel('Indegree d_{in}^i');
    matlab2tikz('figures/mean_indegree_zipf_T=inf.tikz');

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
    
