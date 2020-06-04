function closure_probability = closure_probability(n)

closure_probability = zeros(n,n);

for i=1:n
    for j=1:i-1
        closure_probability(i,j) = 1/i + 1/j - 1/(i*j);
    end
end

closure_probability = closure_probability+closure_probability';

