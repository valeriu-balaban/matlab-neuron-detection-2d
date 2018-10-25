load '~/Downloads/neurons-step1-result-9.mat'
P = cell2mat(NeuronLocations{11});

% Change coordinates from x, y to row, col
TMP = P(:, 1);
P(:, 1) = P(:, 2);
P(:, 2) = TMP;

I = imread('~/Downloads/neurons-step1-result.tif', 11);

[A, PATHS] = TraceConnections(I, P, 15);

%% Add circles

h = viscircles([P(:, 2), P(:, 1)], P(:, 3));

%% Add numbers

for i = 1:size(P, 1)
    text(P(i, 2), P(i, 1), sprintf('%d', i), 'FontSize', 14, 'Color', 'red')
end