load 'three_image.mat'
IMG_STR = 'neurons-step1-result.tif';

num_image = 3;
A = cell(num_image, 1);
PATHS = A;

for image_index = 1:num_image
    
    P = cell2mat(NeuronLocations{image_index});
    % Change coordinates from x, y to row, col
    TMP = P(:, 1);
    P(:, 1) = P(:, 2);
    P(:, 2) = TMP;

    I = imread(IMG_STR, image_index * 10 - 9);

    [A{image_index}, PATHS{image_index}] = TraceConnections(I, P, 15);

    %% Add circles

    h = viscircles([P(:, 2), P(:, 1)], P(:, 3));

    %% Add numbers

    for i = 1:size(P, 1)
        text(P(i, 2), P(i, 1), sprintf('%d', i), 'FontSize', 14, 'Color', 'red')
    end
end
