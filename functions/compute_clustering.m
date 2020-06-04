function clustering = compute_clustering(A)
    n = size(A,2);
    clustering = zeros(n,1);
    triangles  = zeros(n,1);
    closed_triangles = zeros(n,1);
    for i=1:n
        for j=1:n
            if (A(i,j) == 1)
                for k=j+1:n
                    if (A(i, k) == 1)
                        triangles(i) = triangles(i)+1;
                        if (A(j,k) == 1 || A(k,j) == 1)
                            closed_triangles(i) = closed_triangles(i)+1;
                        end
                    end
                end
            end
        end
        if (triangles(i)>0)
            clustering(i) = closed_triangles(i)/triangles(i);
        end
    end  
end