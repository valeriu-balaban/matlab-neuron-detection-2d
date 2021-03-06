function [connected, synapse] = Astar_Synapse(BW, Source, R, Neurons, thetas, ...
                                    fill_gap, theta_thre)                              
% exclude itself
circle_area = circle_points(Neurons, round(R*1.3), BW);
Neurons(Source, :) = [];
R(Source, :) = [];

num = size(Neurons, 1);
% radius for arrival
R_compare = 1.2 * R;
% init output
connected = zeros(num, 1);
synapse = cell(num, 1);
reach_points = zeros(num, 2);
%% Init the Closed list

CLOSED = zeros(0, 2);
CLOSED_COUNT = size(CLOSED, 1);
%% put starting nodes in OPEN list
%OPEN LIST FORMAT
%--------------------------------------------------------------------------
%IS ON LIST 1/0 |X val |Y val |Parent index |h(n) |g(n)|f(n)| target_index
%former theta |
% get start points with source neuron and radius
start_points = circle_area{Source};
index = sub2ind(size(BW), start_points(:,2), start_points(:,1));
start_points = start_points(BW(index), :);
start_points = merge_random(start_points, 10);

% path cost = 0
path_cost = 0;
% put starting nodes in OPEN list
OPEN_COUNT = size(start_points, 1);
OPEN = zeros(OPEN_COUNT, 9);
for i = 1 : OPEN_COUNT
    % parent_index = 0 stands for start
    parent_index = 0;
    % distance
    [goal_distance, target_index] = distance(start_points(i, :), Neurons);
    % theta
    former_theta = thetas(start_points(i, 2), start_points(i, 1));
    % insert
    OPEN(i, :) = insert_open(start_points(i, 1), start_points(i, 2), ...
                             parent_index, path_cost, goal_distance, ...
                             goal_distance, target_index, former_theta);
end
%% Loop 
while true
%% Find out the node with the smallest fn 
    index_min_node = min_fn(OPEN);
    if (index_min_node ~= -1)
    %Set xNode and yNode to the node with minimum fn
        current_Node = OPEN(index_min_node, :);
        current = current_Node(2:3);
        path_cost = current_Node(5);%Update the cost of reaching the parent node
        %Move the Node to list CLOSED
        parent_index = index_min_node;
        CLOSED_COUNT = CLOSED_COUNT + 1;
        CLOSED(CLOSED_COUNT, :) = current;
        OPEN(index_min_node, 1) = 0;
%        plot(current(1), current(2), 'g+');
    else
        %No path exists to the Target!!
        break;
    end%End of index_min_node check
%% whether reach other neurons
    [dis, neu_index] = distance(current, Neurons);
    if dis < R_compare(neu_index)
        real_index = neu_index + (neu_index >= Source);
        if ~connected(neu_index)
            reach_points(neu_index, :) = current;
            connected(neu_index) = current_Node(7);
            fprintf('from neuron %d reach neuron %d\n',Source, real_index);
        elseif current_Node(7) < connected(neu_index)
            reach_points(neu_index, :) = current;
            connected(neu_index) = current_Node(7);
            fprintf('from neuron %d update path to neuron %d\n',Source, real_index);
        end
    else
        former_theta = current_Node(9);
        exp_array = expand_array(current, path_cost, Neurons, CLOSED,...
                           BW, former_theta, thetas, fill_gap, theta_thre);
    end
    exp_count = size(exp_array,1);
    %UPDATE LIST OPEN WITH THE SUCCESSOR NODES
    %OPEN LIST FORMAT
    %--------------------------------------------------------------------------
    %IS ON LIST 1/0 |X val |Y val |Parent index |h(n) |g(n)|f(n)| target_index
    %former theta |
    %--------------------------------------------------------------------------
    %EXPANDED ARRAY FORMAT
    %--------------------------------
    %| X val | Y val | h(n) | g(n) | f(n) | target_index | former_theta
    %--------------------------------
    for i = 1 : exp_count
        % whether expanded point in OPEN list
        exp_point = exp_array(i, :);
        [Isopen, Location] = ismember(exp_point(1:2), OPEN(:,2:3), 'rows');
        if Isopen
            if exp_point(5) < OPEN(Location, 7)
                OPEN(Location, 4:end) = [parent_index, exp_point(3:end)];
            end
        else
            OPEN_COUNT = OPEN_COUNT + 1;
            OPEN(OPEN_COUNT, :) = [1, exp_point(1:2), parent_index, ...
                                   exp_point(3:end)];
        end%End of insert new element into the OPEN list
    end%End of i for
end%End of While Loop

%Once algorithm has run The optimal path is generated by starting of at the
%last node(if it is the target node) and then identifying its parent node
%until it reaches the start node.This is the optimal path

if ~sum(connected)
    fprintf('Fail to find path from neuron %d\n', Source);
else
    reach_index = find(connected);
    for i = 1:sum(connected > 0)
        index = reach_index(i);
        % find last point in each path
        terminal = reach_points(index, :);
        optimal_path = terminal;
        path_len = 1;
        node_index = ismember(OPEN(:, 2:3), terminal, 'rows');
        % parent node of terminal
        parent = OPEN(node_index, 4);
        % path
        while parent ~= 0
            % add current point
            path_len = path_len + 1;
            optimal_path(path_len, :) = OPEN(parent, 2:3);
            % find parent of current point
            parent = OPEN(parent, 4);
        end
        synapse{index} = optimal_path;
    end
end

%% restore indexes including source

connected = [connected(1:Source-1); 0; connected(Source:end)];
synapse = [synapse(1:Source-1); cell(1); synapse(Source:end)];