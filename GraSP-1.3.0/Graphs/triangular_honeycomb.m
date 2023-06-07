%this function gnerates a triangular honeycomb graph with edge length 2
function graph = triangular_honeycomb(N)
    graph = grasp_struct;
    graph.A = zeros(total_nodes(N));
    graph.layout = zeros(total_nodes(N), 2);
    graph.A(1, 2:7) = 1;
    graph.A(total_nodes(N), total_nodes(N - 1) + 1) = 1;
    graph.layout(1,:) = [0 0];
    for layer = 1:N
        next_layer_node = total_nodes(layer) + 1;
        hop = 0;
        x_position = layer;
        y_position = layer*sqrt(3);
        position_add = [1 -sqrt(3); -1 -sqrt(3); -2 0; -1 sqrt(3); 1 sqrt(3); 2 0];
        if layer < N
            graph.A(total_nodes(layer - 1) + 1, total_nodes(layer + 1)) = 1;
            for n = total_nodes(layer - 1) + 1 : total_nodes(layer)
                if n < total_nodes(layer)
                    graph.A(n, n + 1) = 1;
                end
                if ismember(n, total_nodes(layer - 1) + 1 + layer.*(0:5))
                    graph.A(n, next_layer_node - 1 : next_layer_node + 1) = 1;
                    next_layer_node = next_layer_node + 2;
                    hop = hop + 1;
                else
                    graph.A(n, next_layer_node - 1 : next_layer_node) = 1;
                    next_layer_node = next_layer_node + 1;
                end
                graph.layout(n, 1) = x_position;
                graph.layout(n, 2) = y_position;
                x_position = x_position + position_add(hop, 1);
                y_position = y_position + position_add(hop, 2);
            end
        else
            for n = total_nodes(layer - 1) + 1 : total_nodes(layer)
                if n < total_nodes(layer)
                    graph.A(n, n + 1) = 1;
                end
                if ismember(n, total_nodes(layer - 1) + 1 + layer.*(0:5))
                    %next_layer_node = next_layer_node + 1;
                    hop = hop + 1;
                else
                    %next_layer_node = next_layer_node + 2;
                end
                graph.layout(n, 1) = x_position;
                graph.layout(n, 2) = y_position;
                x_position = x_position + position_add(hop, 1);
                y_position = y_position + position_add(hop, 2);
            end
        end
    end
    graph = grasp_symetrise_unweighted(graph);
end

function s = single_layer_nodes(n)
    s = 6*n;
end

function t = total_nodes(layer)
    t = 1+layer*(layer+1)/2*6;
end
