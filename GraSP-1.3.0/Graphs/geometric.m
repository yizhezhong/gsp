%this function gnerates a random geometric graph with N nodes, diameter D, radius r
function graph = geometric(N, D, r)
    graph = grasp_struct;
    graph.A = zeros(N);
    %graph.layout = zeros(total_nodes(N), 2);
    graph.layout = -D + (D+D)*rand(N, 2);

    for i = 1:N
        for j = i:N
            if (graph.layout(i,1) - graph.layout(j,1))^2 + (graph.layout(i,2) - graph.layout(j,2))^2 <= r^2
                graph.A(i,j) = 1;
            end
        end
    end
    graph = grasp_symetrise_unweighted(graph);
end