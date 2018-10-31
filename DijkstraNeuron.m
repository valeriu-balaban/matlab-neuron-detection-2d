function [A, AE, PATHS] = DijkstraNeuron(I, E, U, V, P, N, Q, MAX_FILLS)
%DIJKSTRANEURON Summary of this function goes here
%   Detailed explanation goes here

A       = inf(size(P, 1), 'single');
AE      = inf(size(P, 1), 'single');
PATHS   = cell(size(P, 1));

[m, n]  = size(I);

%% Parameters

GAMMA = 0.8;        % weight of the eigenvalues

THETA = 0.2;

%% Dijkstra Shortest Path

% relative position of the neighbor cells      
delta = [ 0, -1, -1, -1,  0,  1,  1,  1;...
         -1, -1,  0,  1,  1,  1,  0, -1]';

DEST   = cell(size(P, 1));

for i = N
    
    % Distance from source
    C      = inf(m, n, 'single');
    CE     = inf(m, n, 'single'); % numbe of empty pixels in the path

    % Column and Row coordinate of the previous node on the optimal path
    C_prev = zeros(m, n, 'uint16');
    R_prev = zeros(m, n, 'uint16'); 

    % Source
    s_c = P(i, 2);
    s_r = P(i, 1);

    % row index, column index, cost, number of empty pixels in the path
    QUEUE = [s_r, s_c, 0, 0];

    % Initialize
    C(s_r, s_c)  = 0;
    CE(s_r, s_c) = 0;

    % queue is not empty
    while size(QUEUE, 1) > 0
        
        inside = false;

        % find the minimum number of empty pixels in a path
        v   = min(QUEUE(:, 4));
        sel = QUEUE(:, 4) == v;

        % find the minimum cost with the number of empty pixels in a path
        [~, idx] = min(QUEUE(sel, 3));

        % transform the relative index into absolute
        indices = find(sel);
        idx = indices(idx);

        % stop if all of the possible paths have larger gaps
        if v > MAX_FILLS
            send(Q, i);
            break;
        end
        p = QUEUE(idx, 1:2);

        QUEUE(idx, :) = [];        
        
        idx = find(abs(P(:, 1) - p(1)) <= P(:, 3)  &  abs(P(:, 2) - p(2)) <= P(:, 3));
        
        % inside a neuron and it's not the starting one?
        if numel(idx) > 0 && ~ismember(i, idx)
            
            % Update the weight
            for j = idx
                if AE(i, j) >= CE(p(1), p(2)) && A(i, j) > C(p(1), p(2))
                    
                    AE(i, j) = CE(p(1), p(2));
                    %AE(j, i) = CE(p(1), p(2));
                    A(i, j) = C(p(1), p(2));
                    %A(j, i) = C(p(1), p(2));
                    
                    DEST{i, j} = p;
                    %DEST{j, i} = p;
                end
            end
            
            continue;
        end

        if ismember(i, idx)
            inside = true;
        end
        
        % for all neighbors
        for k = 1:size(delta, 1)

            r = p + delta(k, :);

            % Skip if outside image
            if r(1) < 1 || r(1) > m || r(2) < 1 || r(2) > n
                continue;
            end
            
            vx = [U(p(1), p(2)); V(p(1), p(2))];
            vy = [U(r(1), r(2)); V(r(1), r(2))];

            W1 = + GAMMA;
            W2 = + (1 - GAMMA) * norm(delta(k, :)) / 2;

            % compute the cost
            d = delta(k, :) / norm(delta(k, :));
            cost = W1 * (1 - E(p(1), p(2))) + ...
                   W2 * (sqrt(1 - abs(d * vx)) + sqrt(1 - abs(d * vy))) + ...
                   norm(delta(k, :)) - 1;

            num_empty = CE(p(1), p(2)) + (I(r(1), r(2)) < THETA);
            
            if inside
                cost = norm(delta(k, :)) - 1 + (I(r(1), r(2)) < THETA);
                num_empty = 0;
            end

            % update only if less of equal number of empty pixels in the path
            if CE(r(1), r(2)) > num_empty

                CE(r(1), r(2))     = num_empty;
                C(r(1), r(2))      = C(p(1), p(2)) + cost * (I(r(1), r(2)) > THETA);
                R_prev(r(1), r(2)) = p(1);
                C_prev(r(1), r(2)) = p(2);

                % Add it to the list maybe we can decrese the total path
                QUEUE = [QUEUE; r(1), r(2), C(r(1), r(2)), CE(r(1), r(2))];

            elseif CE(r(1), r(2)) == num_empty

                % select the lowest cost
                if C(r(1), r(2)) > cost + C(p(1), p(2))
                    CE(r(1), r(2))     = num_empty;
                    C(r(1), r(2))      = C(p(1), p(2)) + cost * (I(r(1), r(2)) > THETA || inside);
                    R_prev(r(1), r(2)) = p(1);
                    C_prev(r(1), r(2)) = p(2);

                    % Add it to the list maybe we can decrese the total path
                    QUEUE = [QUEUE; r(1), r(2), C(r(1), r(2)), CE(r(1), r(2))];
                end
            end
        end
    end
    
    
    for j = 1:size(P, 1)
        
        if ~isempty(DEST{i, j})
            path = [];
            
            r = DEST{i, j};
            
            % Obtain paths
            while ~all(P(i, 1:2) == r)
                path = [path; r];

                % Move back
                p    = r;
                r(1) = R_prev(p(1), p(2));
                r(2) = C_prev(p(1), p(2)); 
            end
            
            PATHS{i, j} = path;
            %PATHS{j, 1} = path;
        end
    end
end

end

