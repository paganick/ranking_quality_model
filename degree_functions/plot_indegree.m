function plot_indegree(indegree)
    figure();
    hold on;
    for i=1:size(indegree,1)
        plot(indegree(i,:));
    end
    xlabel('t');
    ylabel('Indegree');
    title('Indegree time evolution');
    set(gca, 'XScale', 'log');