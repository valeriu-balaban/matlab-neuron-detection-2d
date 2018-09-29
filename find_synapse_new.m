function [connected, synapse] = find_synapse_new(origins, Neurons, R, source,...
                            skeleton, thetas, theta_thre, fill_gap)
% current: current pixel on the path (start point)
% former: former pixel on the path (tell the forward direction)
% source: the neuron(index) this path originates from 
% Neurons & R: decide the terminal of path(get other neurons / dead path)

%% Variable Definition
% num of Neurons(R)
num = length(R);

% init target neuron
connected = zeros(num, 1);
synapse = cell(num, 1);

num_bifur = size(origins, 1);
% parameters for each bifur 
% current point, former point, path length, path
bifur_paras = cell(num_bifur, 1);
CLOSE_list = cell(num_bifur, 1);

for i = 1 : num_bifur
    bifur_paras{i} = zeros(4, 2);
    bifur_paras{i}(3, 2) = 1;
    bifur_paras{i}(1, :) = origins(i, :); 
    bifur_paras{i}(2, :) = Neurons(source, :);
    bifur_paras{i}(4, :) = Neurons(source, :);
end

[M, N] = size(skeleton);

% figure;
% imshow(skeleton); hold on;

% for compare
% if line is going through the source neuron, stop it
% to prevent mis-stop, set R(source) smaller
R_compare = 1.25 * R;
R_compare(source) = 1.15 * R(source);

while num_bifur > 0
    % get information from bifur parameters
    current = bifur_paras{1}(1, :);
    former = bifur_paras{1}(2, :);
    path = bifur_paras{1}(4:end, :);
    path_length = bifur_paras{1}(3, 1);
    current_CLOSE = CLOSE_list{1};
    % advance along path 1
    while true
        % cover passed points 
        current_CLOSE = [current_CLOSE; current];
        path = [path; current];
        % plot synapses
        % plot([current(1), former(1)],[current(2), former(2)],'LineWidth',2,'Color','blue');
        %% whether reach other neurons
        dis = Neurons - repmat(current, num, 1);
        dis = sqrt(sum(dis.^2, 2));
        if max(dis < R_compare)
            target = find(dis < R_compare);
            num_bifur = num_bifur - 1;
            bifur_paras(1) = []; 
            CLOSE_list(1) = [];
            if target == source
                break;
            end
            path = [path; Neurons(target, :)];
            if connected(target) == 0
                connected(target) = path_length;
                synapse{target} = path;
            elseif path_length < connected(target)
                connected(target) = path_length;
                synapse{target} = path;
            end
            break;
        end 
        % direction of last movement
        former_theta = atan2d(current(2) - former(2), current(1) - former(1));
        %% whether there are bifurcations
        % get neighbours
        [neigh_x, neigh_y] = produce_meshgrid(current, 1, M, N);
        % neighbours which are bright
        index = sub2ind([M, N], neigh_y, neigh_x);
        neigh_index = find(skeleton(index)); 
        neigh_x = neigh_x(neigh_index);
        neigh_y = neigh_y(neigh_index);
        % neighbours which are not in CLOSE
        num_neigh = numel(neigh_x);
        neigh_x = reshape(neigh_x, [num_neigh 1]);
        neigh_y = reshape(neigh_y, [num_neigh 1]);
        
        IsClose = ismember([neigh_x, neigh_y], current_CLOSE, 'rows');
        neigh_x(IsClose) = [];
        neigh_y(IsClose) = [];
        % number of bifurcations
        current_bifur = length(neigh_x);
        %% no way to go
        if ~current_bifur
            % just for test
            if ~fill_gap
                num_bifur = num_bifur - 1;
                bifur_paras(1) = [];
                CLOSE_list(1) = [];
                break;
            end
            % whether this is just a small gap ?
            [area_x, area_y] = produce_meshgrid(current, fill_gap, M, N);
            % bright neighbours
            index = sub2ind([M, N], area_y, area_x);
            bright_index = find(skeleton(index));
            area_x = area_x(bright_index);
            area_y = area_y(bright_index);
            % theta from center to neighbours
            num_area = numel(area_x);
            area_x = reshape(area_x, [num_area 1]);
            area_y = reshape(area_y, [num_area 1]);
            area_theta = atan2d(area_y - current(2), area_x - current(1));
            % get angle from neurite detector
            angle = thetas(current(2), current(1));
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
            in_close = ismember([area_x, area_y], current_CLOSE, 'rows');
            % whether thetas of points in the area are suitable
            index = sub2ind([M, N], area_y, area_x);
            neurite_theta = thetas(index);
            theta_suit = abs(neurite_theta - forward_theta);
            theta_suit(theta_suit > 180) = 360 - theta_suit(theta_suit > 180);
            theta_suit = theta_suit < theta_thre;
            % available point
            available = in_direction & ~in_close & theta_suit;
            if sum(available) > 0
                % continue with the point after gap
                dist = sqrt((area_x - current(1)).^2 + (area_y - current(2)).^2);
                dist_available = dist ./ available;
                [dis_min, index_forward] = min(dist_available);
                next_after_gap = [area_x(index_forward), area_y(index_forward)];
                former = current;
                current = next_after_gap;
                path_length = path_length + dis_min;
            else
                % end this path
                num_bifur = num_bifur - 1;
                bifur_paras(1) = [];
                CLOSE_list(1) = [];
                break;
            end
        %% more way to go
        elseif current_bifur > 1
            theta_neigh = atan2d(neigh_y - current(2), neigh_x - current(1));
            % sort directions according to the difference between former
            % direction and forward direction.
            del_theta = abs(theta_neigh - former_theta);
            del_theta(del_theta > 180) = 360 - del_theta(del_theta > 180);
            [~, del_index] = sort(del_theta);
            % num of all bifurs
            num_bifur = num_bifur + current_bifur - 1;
            new_bifur_paras = cell(current_bifur, 1);
            new_CLOSE_list = cell(current_bifur, 1);
            current_CLOSE = [current_CLOSE; neigh_x, neigh_y];
            % value of bifurs
            for i = 1:current_bifur
                new_start = [neigh_x(del_index(i)), neigh_y(del_index(i))];
                step = sqrt(sum((new_start - current) .^ 2));
                new_bifur_paras{del_index(i)} = [new_start; current; ...
                                      path_length + step, 0; [path; new_start]];
                new_CLOSE_list{i} = current_CLOSE;
            end
            bifur_paras(1) = [];
            CLOSE_list(1) = [];
            bifur_paras = [new_bifur_paras; bifur_paras];
            CLOSE_list = [new_CLOSE_list; CLOSE_list];
            break;
        %% one way to go
        else
            former = current;
            current = [neigh_x, neigh_y];
            path_length = path_length + sqrt(sum((current - former) .^ 2));
        end
    end
end
end