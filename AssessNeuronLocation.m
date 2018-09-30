function p = AssessNeuronLocation(I, R_MIN, R_MAX)
% ASSESSNEURONLOCATION Predicts the most probable location and radius for a
% neuron in the gray image I
%
% Input parameters:
%   I               Gray image containing the neuron
%   RADIUS_LIST     List of integers with the radii to be assessed
%
% Return values:
%   p               An 1-by-3 vector with the coordinates of the neuron and
%                   its radius. [COLUMN, ROW, RADIUS]


%% Parameters
% Fixed parameters that depend on the properties of the images containing
% the neurons for detection

% Preprocessing
BW_THRESHOLD        = 45;
INTENSITY_THRESHOLD = 190;
AREA_THRESHOLD      = 500;
EXTENT_THRESHOLD    = 0.4;

EDGE_WIDTH = 3;

%% Image preprocessing

% Step 1.
% -----------------
% Sharpen the image by assigning to each pixel the value of the local 
% max or min, whichever is closer.

H = strel('disk', 2).Neighborhood;

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


I = imclose(I > BW_THRESHOLD, strel('disk', 2));
I = bwareaopen(I, 50);

% Step 4.
% -----------------
% fill holes in neurons which caused when removing dots

I_f = imfill(I, 'holes') & ~I;

% preventing from fill very large holes
I = I | I_f & ~bwareaopen(I_f, 120);


%% Neuron Detection
% For each radius in the list, compute the ratio of pixels at the edge 
% forming different angles with the center. Select the radius that give the
% maximum ratio as the radius of the neuron

RADIUS_LIST = R_MIN:R_MAX;

row = round(size(I, 1) / 2);
col = round(size(I, 1) / 2);

rows = max(1, row - R_MAX - EDGE_WIDTH) : min(size(I, 1), row + R_MAX + EDGE_WIDTH);
cols = max(1, col - R_MAX - EDGE_WIDTH) : min(size(I, 2), col + R_MAX + EDGE_WIDTH);

Img = I(rows, cols);

theta  = 0;
radius = R_MIN;

% Precompute the distance and angle matrix when the dimensions changed
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

for i = 1:numel(RADIUS_LIST)

    % Compute the ratio of angles found on the edge compared to number
    % of posible angles for that radius
    edge_angles = zeros(1, 361, 'uint8');
    edge_angles(angles(edge_pixels(:, :, i) & Img) + 181) = 1;
    ratio       = sum(edge_angles) / num_angles(i);

    % Save the radius with the maximum ratio
    if ratio > theta
        theta  = ratio;
        radius = RADIUS_LIST(i);
    end
end

p = [col, row, radius];

end

