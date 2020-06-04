function p = indegree_prob(N,d,i,T)
p = 0;
p2 = 1/i*(1-((N-i-1)/(N-1))^T);
if (i==1)
    b = nchoosek(N-1, d);
    p = b*p2^d*(1-p2)^(N-i-d);
else
    for k=0:d
        if (k>i-1)
            b1 = 0;
        else
            b1 = nchoosek(i-1,k);
        end
        if (d-k>N-i)
            b2 = 0;
        else
            b2 = nchoosek(N-i, d-k);
        end
        p1 = 1/(i-1)*(1-((N-i)/(N-1))^T);
        p = p + b1*p1^k*(1-p1)^(i-1-k)*b2*p2^(d-k)*(1-p2)^(N-i-d+k);
    end
end
end
