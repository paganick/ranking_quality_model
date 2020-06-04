function p = outdegree_prob(N) % can be optimized!
p=zeros(N-1,N-1);
p(1,:)=1/(N-1)*ones(1,N-1);
for d=2:N-1
    for k_0=1:N-1
        for k=k_0:N-d
            p(d,k_0) = p(d,k_0)+1/(N-1-k)*p(d-1,k+1);
        end
    end
end
end