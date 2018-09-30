function exp_array = expand_array(Node, hn, Neurons, CLOSED, BW, former_theta,...
                                  thetas, fill_gap, theta_thre)
    %Function to return an expanded array
    %This function takes a node and returns the expanded list
    %of successors,with the calculated fn values.
    %The criteria being none of the successors are on the CLOSED list.
    
    [M, N] = size(BW);
    % meshgrid
    [x_grid, y_grid] = produce_meshgrid(Node, 1, N, M);
    % get white points
    index = sub2ind([M, N], y_grid, x_grid);
    bright_index = find(BW(index));
    x_grid = x_grid(bright_index);
    y_grid = y_grid(bright_index);
    
    % whether neighbour points in closed list
    exp_points = [x_grid, y_grid];
    Isclose = ismember(exp_points, CLOSED, 'rows');
    exp_points(Isclose, :) = [];
    
    exp_array = exp_points;
    exp_count = size(exp_array, 1);
    
    if ~exp_count && fill_gap > 0
        % whether this is just a small gap ?
        [area_x, area_y] = produce_meshgrid(Node, fill_gap, M, N);
        % bright neighbours
        index = sub2ind([M, N], area_y, area_x);
        bright_index = find(BW(index));
        area_x = area_x(bright_index);
        area_y = area_y(bright_index);
        % theta from center to neighbours
        area_theta = atan2d(area_y - Node(2), area_x - Node(1));
        % get angle from neurite detector
        angle = thetas(Node(2), Node(1));
        % since angle means two opposite direction, 
        % we should find the direction to forward
        theta1 = angle; 
        theta2 = - sign(theta1) * (180 - abs(angle));
        theta_dis = abs([theta1, theta2] - former_theta);
        theta_dis(theta_dis > 180) = 360 - theta_dis(theta_dis > 180);
        [~, index_theta] = min(theta_dis);
        forward_theta = (index_theta == 1) * theta1 + (index_theta == 2) * theta2;
        % average of neurite detector and last movement
        forward_theta = mean([forward_theta, former_theta]); 
        theta_differ = abs(area_theta - forward_theta);
        % whether points of the area are in forward direction
        in_direction = theta_differ < theta_thre;
        % whether points of the area are in CLOSE
        in_close = ismember([area_x, area_y], CLOSED, 'rows');
        % whether thetas of points in the area are suitable
        index = sub2ind([M, N], area_y, area_x);
        neurite_theta = thetas(index);
        theta_suit = abs(neurite_theta - forward_theta);
        theta_suit(theta_suit > 180) = 360 - theta_suit(theta_suit > 180);
        theta_suit = theta_suit < theta_thre;
        % available point
        available = in_direction & ~in_close;
        % available = in_direction & ~in_close & theta_suit;
        if sum(available) > 0
            % continue with the point after gap
            dist = sqrt((area_x - Node(1)).^2 + (area_y - Node(2)).^2);
            dist_available = dist ./ available;
            [~, index_forward] = min(dist_available);
            next_after_gap = [area_x(index_forward), area_y(index_forward)];
            exp_array = next_after_gap;
            exp_count = 1;
        end
    end
    
    % calculate hn, gn and fn
    hn_new = hn + sqrt(sum((repmat(Node, exp_count, 1) - exp_array).^2, 2));
    
    gn_new = zeros(exp_count, 1);
    
    target_index = gn_new;
    
    for i = 1:exp_count
        [gn_new(i), target_index(i)] = distance(exp_array(i, :), Neurons);
    end
        
    fn_new = hn_new + gn_new;
    
    former_theta = atan2d(exp_array(:, 2) - Node(2), exp_array(:, 1) - Node(1));
    
    exp_array = [exp_array, hn_new, gn_new, fn_new, target_index, former_theta];
end  