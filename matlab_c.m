load('~/Downloads/neurons-yunpeng.mat')
IMG_STR = '~/Downloads/neurons-yunpeng.tif';

image_index = 1;

P = cell2mat(NeuronLocations{image_index});

% Change coordinates from x, y to row, col
TMP = P(:, 1);
P(:, 1) = P(:, 2);
P(:, 2) = TMP;

I = imread(IMG_STR, image_index);

[A, PATHS] = TraceConnections(I, P, 25);

%% Add circles

h = viscircles([P(:, 2), P(:, 1)], P(:, 3));

%% Add numbers

for i = 1:size(P, 1)
    text(P(i, 2), P(i, 1), sprintf('%d', i), 'FontSize', 14, 'Color', 'red')
end
