function [tree, n_nodes, expected_clustering_coefficient, cumulative_probability] = add_child(tree, parent_id, max_level, n_nodes, y, expected_clustering_coefficient, cumulative_probability, max_list_size, connection_probability_matrix, min_p)
    level = tree(parent_id).level+1;
    if (level <= max_level)  
        child_id = n_nodes+1;
   %     tree(parent_id).children = [tree(parent_id).children, child_id];
        tree(child_id).level  = level;
        tree(child_id).expected_closed_triangles  = tree(parent_id).expected_closed_triangles;
        if (y == 1)
            tree(child_id).p      = tree(parent_id).p * (1/level);
            for i=1:tree(parent_id).list_size
                element = tree(parent_id).list(i);
                tree(child_id).expected_closed_triangles = tree(child_id).expected_closed_triangles + connection_probability_matrix(element, level);
            end
            tree(child_id).list   = [tree(parent_id).list, level];
            tree(child_id).list_size = tree(parent_id).list_size+1;
            tree(child_id).n_triangles = tree(parent_id).n_triangles + (tree(parent_id).list_size);

        elseif (y==0)
            tree(child_id).p      = tree(parent_id).p * (1-1/level);
            tree(child_id).list   = [tree(parent_id).list];
            tree(child_id).list_size = tree(parent_id).list_size;
            tree(child_id).n_triangles = tree(parent_id).n_triangles;
        end
   %     tree(child_id).n_triangles = (tree(child_id).list_size)*(tree(child_id).list_size-1)/2;
   %     tree(child_id).y      = y;
   %     tree(child_id).parent = parent_id;
   %     if (tree(child_id).n_triangles>0)
        if (tree(child_id).n_triangles>0)
   %         tree(child_id).expected_clustering_coefficient = tree(child_id).expected_closed_triangles/tree(child_id).n_triangles * tree(child_id).p;   
            child_expected_clustering_coefficient = tree(child_id).expected_closed_triangles/tree(child_id).n_triangles * tree(child_id).p;   
            %expected_clustering_coefficient(level,1) = expected_clustering_coefficient(level,1) + tree(child_id).n_triangles* tree(child_id).p;
            %expected_clustering_coefficient(level,2) = expected_clustering_coefficient(level,2) + tree(child_id).expected_closed_triangles* tree(child_id).p;
            expected_clustering_coefficient(level,3) = expected_clustering_coefficient(level,3) + child_expected_clustering_coefficient;
        end
        n_nodes = child_id;
        cumulative_probability(level) = cumulative_probability(level) + tree(child_id).p;
        [tree, ~, expected_clustering_coefficient, cumulative_probability] = add_child(tree, child_id, max_level, n_nodes, 0, expected_clustering_coefficient, cumulative_probability, max_list_size, connection_probability_matrix, min_p);
        if (tree(child_id).list_size<=max_list_size)
            if (tree(child_id).p>min_p)
                [tree, ~, expected_clustering_coefficient, cumulative_probability] = add_child(tree, child_id, max_level, n_nodes, 1, expected_clustering_coefficient, cumulative_probability, max_list_size, connection_probability_matrix, min_p);
            end
        end
    end
end
