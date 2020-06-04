function cost = cost(A,i)
    if (sum(A(i,:)) > 0)
%        cost = sum(A(i,:))*log2(sum(A(i,:)));
        cost = sum(A(i,:));
    else
        cost = 0;
    end
end

