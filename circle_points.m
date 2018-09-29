function circle_area = circle_points(neurons, R, I)

if size(neurons, 1) == 0
    return
end

[M, N] = size(I);
num = length(R);
circle_area = cell(num, 1);

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
    circle_area{k} = [area_x(in_circle), area_y(in_circle)];
end