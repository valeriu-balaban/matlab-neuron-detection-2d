function [position, radius] = FindNeurons(I, RADIUS_LIST, EDGE_WIDTH, THETA_THRESHOLD)
% FINDNEURONS Implements image processing techniques for extracting
% neuron locations and sizes from gray images
%
% Input parameters:
%   I               Gray image containing the neurons
%   RADIUS_LIST     List of integers indicating the radii of the neurons to
%                   be checked for detectio
%   EDGE_WIDTH      The width in pixels of the neuron edge used to compute
%                   the THETA_THRESHOLD
%   THETA_THRESHOLD The arc ratio defined by the edge pixels. Set to 1 for
%                   360 degrees
%
% Return values:
%   position        An N-by-2 matrix with the coordinates of the neuron
%                   centers found
%   radius          An N-vector containing the estimated radius in pixels

%% Parameters
% Fixed parameters that depend on the properties of the images containing
% the neurons for detection

% Preprocessing
BW_THRESHOLD        = 45;
INTENSITY_THRESHOLD = 190;
AREA_THRESHOLD      = 500;
EXTENT_THRESHOLD    = 0.4;

% Detection
DENSITY_LIST        = 0.3:0.05:1.0;
R_SOMA_SCALE        = 0.6;
SOMA_THRESHOLD      = 0.4;

%% Image preprocessing

wh = waitbar(0, 'Preprocessing the image...', ...
                'Name', 'Neuron Detection in Progress', ...
                'WindowStyle', 'modal');

% Step 1.
% -----------------
% Sharpen the image by assigning to each pixel the value of the local 
% max or min, whichever is closer.

H = strel('disk', 1).Neighborhood;

I_max = ordfilt2(I, sum(H(:)), H);
I_min = ordfilt2(I, 1, H);

I_idx = I_max - I > I - I_min;

I = I_max;
I(I_idx) = I_min(I_idx);


% Step 2.
% -----------------
% Remove small neurons with no connection that appear on the images as 
% small white dots.

% Label regions so we can filter them later
L = bwlabel(I > BW_THRESHOLD);

stats = regionprops(L, I, 'Extent', 'Area', 'MeanIntensity'); %#ok<MRPBW>
idx   = find([stats.MeanIntensity] > INTENSITY_THRESHOLD & ...
             [stats.Area] < AREA_THRESHOLD & [stats.Extent] > EXTENT_THRESHOLD);

% Set the pixels of the regions that satisfy the above criteria to 0
I(ismember(L, idx)) = 0;


% Step 3.
% -----------------
% Prepare the black and white image for neuron detection.
% Remove small unconnected white patches on the image and slightly extend
% everything else.

BW = I > BW_THRESHOLD;

BW = imclose(BW, strel('disk', 2));

BW = bwareaopen(BW, 50);

% Step 4.
% -----------------
% fill holes in neurons which caused when removing dots

BW_filled = imfill(BW, 'holes');
fill_differ = BW_filled & ~BW;
% preventing from fill very large holes
fill = fill_differ & ~bwareaopen(fill_differ, 120);
BW = BW | fill;

I = I .* uint8(BW);

% Free memory
clearvars  I_max I_min I_idx L;
%% Neuron Detection

% Step 1.
% -----------------
% Convolve the image with disk of different radius that match the neuron
% sizes. This produces peaks at locations where the disk covers the most
% number of white pixels. Use theresholding to find regions containing
% peaks and add their centroid to the list of potential neurons.

% Remark
% -----------------
% Potential improvement is to replace the centroid with the local maximas
% in the image which are above a threshold.

points = cell(1, numel(RADIUS_LIST) * numel(DENSITY_LIST));
k = 1;

for radius = RADIUS_LIST
    waitbar(0.1 + 0.45 * k / (numel(RADIUS_LIST) * numel(DENSITY_LIST)), ...
            wh, sprintf('Finding neurons of radius %i...', radius));
    
    I_d = conv2(BW, single(fspecial('disk', radius)), 'same');

    for T = DENSITY_LIST
        stats = regionprops(I_d > T, 'Centroid');
        
        points{k} = reshape([stats.Centroid], 2, []);
        k = k + 1;
    end
end

points = round(transpose(cell2mat(points)));

% Free memory
clearvars I_d;


% Step 2.
% -----------------
% For each posible point and radius in the list, check if the ratio of 
% pixels at the edge forming a different angle with the center is above the
% specified threshold

% Output
radius = zeros(size(points,1), 1);
theta  = zeros(size(points,1), 1);


% Intermediate varibles
R_MAX       = RADIUS_LIST(end);
D           = [];
num_angles  = zeros(1, numel(RADIUS_LIST));


for k = 1:size(points,1)
    if mod(k, 79) == 0
        waitbar(0.55 + 0.45 * k / size(points,1), wh,...
                sprintf('Evaluating neuron %i of %i...', k, size(points,1)));
    end
    
    row = points(k, 2);
    col = points(k, 1);
    
    rows = max(1, row - R_MAX - EDGE_WIDTH) : min(size(I, 1), row + R_MAX + EDGE_WIDTH);
    cols = max(1, col - R_MAX - EDGE_WIDTH) : min(size(I, 2), col + R_MAX + EDGE_WIDTH);
    
    Img = BW(rows, cols);
    
    % Recompute the distance and angle matrix when the dimensions changed
    if size(D, 1) ~= numel(rows) || size(D, 2) ~= numel(cols)
        [mc, mr]    = meshgrid(cols - col, rows - row);
        D           = sqrt(single(mc .^ 2 + mr .^ 2));
        angles      = int16(atan2d(mc, mr));
        
        edge_pixels = zeros(size(D, 1), size(D, 2), numel(RADIUS_LIST), 'logical');
        
        for i = 1:numel(RADIUS_LIST)
            
            % For each radius compute the maximum number of angles that can
            % be found inside
            inside_angles = zeros(1, 361, 'uint8');
            inside_angles(angles(D <= RADIUS_LIST(i)) + 181) = 1; 
            num_angles(i) = sum(inside_angles);
            
            % Create the mask of edge pixels for each radius
            edge_pixels(:, :, i) = (D <= RADIUS_LIST(i) + EDGE_WIDTH) & ...
                                   (D >= RADIUS_LIST(i) - EDGE_WIDTH);
        end
    end
    
    for i = 1:numel(RADIUS_LIST)
        
        % Compute the ratio of angles found on the edge compared to number
        % of posible angles for that radius
        edge_angles = zeros(1, 361, 'uint8');
        edge_angles(angles(edge_pixels(:, :, i) & Img) + 181) = 1;
        ratio       = sum(edge_angles) / num_angles(i);
        
        % Save the radius with the maximum ratio
        if ratio > theta(k)
            theta(k)  = ratio;
            radius(k) = RADIUS_LIST(i);
        end
    end
    
    % compute the bright ratio in soma area
    radius_soma = radius(k) * R_SOMA_SCALE;
    rows_soma = max(1, row - radius_soma) : min(size(I, 1), row + radius_soma);
    cols_soma = max(1, col - radius_soma) : min(size(I, 2), col + radius_soma);
    [mcs, mrs] = meshgrid(cols_soma, rows_soma);
    % bright num in circle / num in circle
    D = sqrt(single((mcs - col) .^ 2 + (mrs - row) .^ 2));
    in_soma = D < radius_soma;
    img_soma = BW(round(rows_soma), round(cols_soma));
    soma_ratio = sum(sum(img_soma & in_soma)) / sum(in_soma(:));
    
    % if not meet the need of soma, theta = 0
    if soma_ratio < SOMA_THRESHOLD
         theta(k) = 0;
    end
end

% Select only points above the threshold
idx = theta > THETA_THRESHOLD;

position = points(idx, :);
radius  = radius(idx);
theta   = theta(idx);

% Step 3.
% -----------------
% Remove neighbors with smaller theta within R_MAX 

for i = 1:size(position, 1)
    
    if isnan(position(i, 1))
        continue;
    end
    
    d   = sqrt((position(:, 1) - position(i, 1)) .^ 2 + ...
               (position(:, 2) - position(i, 2)) .^ 2);
    idx = (d < R_MAX) & (theta <= theta(i));
    
    % Prevent romoving current point
    idx(i) = 0;

    position(idx,:) = NaN;
end

idx = ~isnan(position(:, 1));

radius   = radius(idx);
position = position(idx, :);

close(wh);
end