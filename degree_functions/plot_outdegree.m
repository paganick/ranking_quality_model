function plot_outdegree(outdegree)
    figure();
    hold on;
    for i=1:size(outdegree,1)
        plot(outdegree(i,:));
    end
    xlabel('t');
    ylabel('Outdegree');
    title('Outdegree time evolution');
    set(gca, 'XScale', 'log');
