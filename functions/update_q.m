function q = update_q(A,q,hub)
    p_bad = 1;
    for j=1:size(A,1)
       if (A(hub, j) == 1)
           p_bad = p_bad*(1-q(j));
       end
    end
    q(hub) = 1 -p_bad;
end

