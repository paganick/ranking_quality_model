function A = SF_simulate(model)
    edges = preferentialAttachment(model.n,model.outdegree);
    A = sparse(edges(:,1),edges(:,2),1, model.n, model.n);
    A = full(A);