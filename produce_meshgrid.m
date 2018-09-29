function [x_grid, y_grid] = produce_meshgrid(center, r, m, n)
    % produce meshgrid with center and range distance
    % x_range
    x = center(1);
    y = center(2);
    x_range = x - r : x + r;
    exceed_x = x_range < 1 | x_range > n;
    x_range(exceed_x) = [];
    % y_range
    y_range = y - r : y + r;
    exceed_y = y_range < 1 | y_range > m;
    y_range(exceed_y) = [];
    [x_grid, y_grid] = meshgrid(x_range, y_range);  
end