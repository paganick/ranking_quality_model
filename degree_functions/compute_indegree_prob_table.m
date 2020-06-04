function indegree_prob_table = compute_indegree_prob_table(model)
    indegree_prob_table = zeros(model.n-1, model.n);
    for d=0:model.n-1
        for k=1:model.n
            indegree_prob_table(k, d+1) = indegree_prob(model.n,d,k,model.T);
        end
    end
end