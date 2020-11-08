function [indegree_prob_table,indegree_prob_table_inf] = compute_indegree_prob_table(model)
    indegree_prob_table     = zeros(model.n, model.n);
    indegree_prob_table_inf = zeros(model.n, model.n);
    nchoosek_table = zeros(model.n, model.n);
    tic
    for i=0:model.n
        for j=0:i
            nchoosek_table(i+1, j+1) = nchoosek(i,j);
        end
    end
    

    for d=0:model.n-1
        for i=1:model.n
            [indegree_prob_table(i, d+1), indegree_prob_table_inf(i, d+1)] = indegree_prob(model.n,d,i,model.T,nchoosek_table);
        end
    end
    toc
end