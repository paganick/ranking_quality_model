function A = BR(A, i, j, q)
    if (A(i,j) == 0)
        u = utility(A,i,q);
        A_alt = A; 
        A_alt(i,j) = 1-A_alt(i,j);   
        if (q(j) >= u)
               A = A_alt;
        end
    end
end

