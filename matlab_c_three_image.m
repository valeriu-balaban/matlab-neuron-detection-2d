load '~/downloads/three_image.mat'
IMG_STR = '~/downloads/three_image.tif';

num_image = 3;
A = cell(num_image, 1);
PATHS = A;

for image_index = 1:num_image
    %% remove previous files
    if exist('data.txt', 'file')
        delete('data.txt');
    end
    if exist('result.txt', 'file')
        delete('result.txt');
    end
    %% compute matrix
    P = cell2mat(NeuronLocations{image_index});
    % Change coordinates from x, y to row, col
    TMP = P(:, 1);
    P(:, 1) = P(:, 2);
    P(:, 2) = TMP;

    I = imread(IMG_STR, image_index);

    [A{image_index}, PATHS{image_index}] = TraceConnections(I, P, 15);

    %% Add circles

    h = viscircles([P(:, 2), P(:, 1)], P(:, 3));

    %% Add numbers

    for i = 1:size(P, 1)
        text(P(i, 2), P(i, 1), sprintf('%d', i), 'FontSize', 14, 'Color', 'red')
    end
end
