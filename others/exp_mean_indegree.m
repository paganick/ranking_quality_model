function d = exp_mean_indegree(model, i)
    d = model.n/i - ( ((model.n-i)/(model.n-1))^model.T + ((model.n-i-1)/(model.n-1))^model.T*(model.n-i)/i);
end
