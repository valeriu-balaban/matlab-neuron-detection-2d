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

% Step 1.
% -----------------
% Convolve the image with disk of different radius that match the neuron
% sizes. This produces peaks at locations where the disk covers the most
% number of white pixels. Select the highest peak as the location and the 
% radius of the neuron



RADIUS_LIST = R_MIN:R_MAX;

value  = zeros(numel(RADIUS_LIST), 1);
center = zeros(numel(RADIUS_LIST), 2);

for k = 1:numel(RADIUS_LIST)
    
    I_d = conv2(I, single(fspecial('disk', RADIUS_LIST(k))), 'same');

    [value(k), idx] = max(I_d(:));
    [center(k, 2), center(k, 1)] = ind2sub(size(I_d), idx);
end

% Use first peak as an indicator of neuron center and size
[~, l] = findpeaks(value .* RADIUS_LIST');

p = [center(l(1), :), RADIUS_LIST(l(1))];

disp(p)

end

