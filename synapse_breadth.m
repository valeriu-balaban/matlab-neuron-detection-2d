function breadth = synapse_breadth(synapse, BW)

num = length(synapse);

breadth = cell(num);


sigma = 4;
[U, V, thetas] = neurite_vector(BW, sigma);
% figure;
% imshow(BW);
% hold on;
% quiver(U, V);

for i = 1 : num-1
    for j = i+1 : num
        if ~isempty(synapse(i, j))
            path = synapse{i, j};
            len = size(path, 1);
            breadth{i, j} = zeros(len, 1);
            for k = 1 : len
                % the range of direction is -180° ~ 180°
                % first perpendicular direction
                perpendicular_theta1 = thetas(path(k, 2), path(k, 1)) + 90;
                % if theta > 180, then theta = theta - 360
                perpendicular_theta1 = perpendicular_theta1 - ...
                                        (perpendicular_theta1 > 180) * 360;
                breadth1 = first_black_in_direction(path(k, :), perpendicular_theta1, BW);
                % second perpendicular direction
                perpendicular_theta2 = thetas(path(k, 2), path(k, 1)) - 90;
                % if theta < -180, then theta = theta + 360
                perpendicular_theta2 = perpendicular_theta2 + ...
                                        (perpendicular_theta2 < -180) * 360;
                breadth2 = first_black_in_direction(path(k, :), perpendicular_theta2, BW);
                breadth{i, j}(k) = breadth1 + breadth2;
            end
        end
    end
end

    function distance = first_black_in_direction(point, theta, BW)
        distance = 1;
        [m, n] = size(BW);
        while true
            point_in_direction = point + distance * [cosd(theta), sind(theta)];
            point_in_direction = round(point_in_direction);
            if (point_in_direction(1) < 1 || point_in_direction(1) > n ||...
               point_in_direction(2) < 1 || point_in_direction(2) > m)
                break
            elseif BW(point_in_direction(2), point_in_direction(1))
                distance = distance + 1;
            else
                break;
            end
        end
        distance = sqrt(sum((point - point_in_direction) .^ 2));
    end
end