function draw_circles(neurons, R, I)

figure;
imshow(I), hold on;
if size(neurons, 1) == 0
    return
end
axis on, xlabel x, ylabel y;
plot(neurons(:,1),neurons(:,2),'.','color','red', 'MarkerSize', 2);

[M, N] = size(I);
num = length(R);

for k = 1:num
    % generate meshgrid
    annulus = 1;
    row = max(1, neurons(k,1) - R(k) - annulus):...
        min(N, neurons(k,1) + R(k) + annulus);
    col = max(1, neurons(k,2) - R(k) - annulus)...
        :min(M, neurons(k,2) + R(k) + annulus);
    [area_x, area_y] = meshgrid(row, col);

    dist = sqrt((area_x - neurons(k,1)).^2 + (area_y - neurons(k,2)).^2);
    % number of points in the circular area of points(k)
    in_circle = dist >= R(k) - annulus & dist <= R(k) + annulus;
    plot(area_x(in_circle), area_y(in_circle), '.','color','red', 'MarkerSize', 1);
end