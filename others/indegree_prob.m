function [p, p_inf] = indegree_prob(N,d,i,T, nchoosek_table)
p2 = 1/i*(1-((N-i-1)/(N-1))^T);
p2_inf = 1/i;
if (i==1)
    b = nchoosek(N-1, d);
    p = b*p2^d*(1-p2)^(N-i-d);
    p_inf = b*p2_inf^d*(1-p2_inf)^(N-i-d);
else
    p = 0;
    p_inf = 0;
    for k=0:d
        if (k>i-1)
            b1 = 0;
        else
            %b1 = nchoosek(i-1,k);
            b1 = nchoosek_table(i, k+1);
        end
        if (d-k>N-i)
            b2 = 0;
        else
            %b2 = nchoosek(N-i, d-k);
            b2 = nchoosek_table(N-i+1, d-k+1);
        end
        p1 = 1/(i-1)*(1-((N-i)/(N-1))^T);
        p1_inf = 1/(i-1);
        if (b1>0 && b2>0)
            p = p + b1*p1^k*(1-p1)^(i-1-k)*b2*p2^(d-k)*(1-p2)^(N-i-d+k);
            p_inf = p_inf + b1*p1_inf^k*(1-p1_inf)^(i-1-k)*b2*p2_inf^(d-k)*(1-p2_inf)^(N-i-d+k);
        end
        %if (p<1.0e-6)
        %    p = 0;
        %end
        %if (p_inf<1.0e-6)
        %    p_inf = 0;
        %end
    end
end
end
