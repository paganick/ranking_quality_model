function A = simulate(model)
    A = zeros(model.n);
    %indegree = zeros(model.n, model.T);
    %outdegree = zeros(model.n, model.T);    
    for t=1:model.T
        old_A = A;
        for i=1:model.n
            j = [];
%            indeg = sum(A,1)' + ones(model.n,1);
            while (not(length(unique(j))==1) || (j==i))
                j = randi(model.n,1,1);
%               j = weightedRandomSample(1,[1:n],indeg/sum(indeg));
            end
            A = BR(A, i, j, model.q);
        end
        %outdegree(:, t) = sum(A,2);
        %indegree(:,t) = sum(A,1)';
        if (A == old_A)
            convergence = is_converged(A, model);
            if (convergence == true)
                break
            end
        end
    end
    %t
    %plot_outdegree(outdegree);
    %plot_indegree(indegree);
    % figure();
    % G = digraph(A);
    % plot(G);
end

function convergence = is_converged(A, model)
    convergence = true;
    if (A(1,2) ~= 1)
        convergence = false;
    else
        for j=2:model.n
            if (A(j,1) ~=1)
                convergence = false;
                break
            end
        end
    end
end
        