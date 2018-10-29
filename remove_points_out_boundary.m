point_1 = [2000, 4000];
point_2 = [1700, 2500];
point_3 = [3000, 1500];
point_4 = [7500, 4000];

% img1 = remove_parts(img, 1, point_1, point_2);
% img2 = remove_parts(img1, 0, point_2, point_3);
% img3 = remove_parts(img2, 0, point_3, point_4);

%% remove points outside the boundaries
load 'neurons-step1-result.mat'
index = 11;
neurons_index = cell2mat(NeuronLocations{index});
neurons_1 = remove_points(neurons_index, 1, point_1, point_2);
neurons_2 = remove_points(neurons_1, 0, point_2, point_3);
neurons_3 = remove_points(neurons_2, 0, point_3, point_4);

figure; imshow(image);
hold on;
h = viscircles([neurons_3(:, 1), neurons_3(:, 2)], neurons_3(:, 3));

NeuronLocations{index} = num2cell(neurons_3);
save('neurons-step1-result.mat', 'NeuronLocations');

%% remove points over the line of p1 and p2
function points = remove_points(points, upordown, p1, p2)
    % upordown: 0 for up and 1 for down
    k = (p1(2) - p2(2)) / (p1(1) - p2(1));
    b = p1(2) - k * p1(1);
    N = size(points, 1);
    stay_flag = true(N, 1);
    for i = 1:N
        if (points(i, 2) * k + b) >= points(i, 1) && ~upordown || ...
           (points(i, 2) * k + b) <= points(i, 1) && upordown
            stay_flag(i) = false;
        end
    end
    points = points(stay_flag, :);
end