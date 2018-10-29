function [A, PATHS] = TraceConnections(I, P, MAX_FILLS)
%TRACECONNECTIONS Summary of this function goes here
%   Detailed explanation goes here


A       = inf(size(P, 1));
PATHS   = cell(size(P, 1));


%% Parameters

SIGMA = 1;          % variance of the convolution kernel
GAMMA = 0.8;        % proportional weight of the eigenvalue when tracing

% Transformation matrix for the eigenvalues
T = [1, 1/3; 1/3, 1];

% Image Preprocessing
BW_THRESHOLD        = 45;
INTENSITY_THRESHOLD = 190;
AREA_THRESHOLD      = 500;
EXTENT_THRESHOLD    = 0.4;


%% Image preprocessing

IMG = im2double(I);

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

%% Hessian computation

I = im2double(I);

gaussian_img = imgaussfilt(I, SIGMA, 'padding', 'replicate');

% second-order derivatives
[gx, gy] = gradient(gaussian_img);
[gxx, gxy] = gradient(gx);
[gyx, gyy] = gradient(gy);


%% Compute line directions
[m, n] = size(I);
E      = zeros(m, n);
U      = zeros(m, n);
V      = zeros(m, n);

for i = 1:m
   parfor j = 1:n
        Hessian = [gxx(i,j) gxy(i,j); gyx(i,j) gyy(i,j)];
        [eigen_vector, eigen_value] = eig(Hessian, 'vector');
        
        % Adjust eigenvalues
        eigen_value = abs(T * eigen_value);
        
        E(i, j) = max(eigen_value);
        [~, idx] = min(eigen_value);
        
        % eigenvector
        U(i, j) = eigen_vector(1, idx);
        V(i, j) = eigen_vector(2, idx);
   end
   
   fprintf(1, 'Finding ridges, line %i of %i\n', i, m);
end

% Normalize eigenvalues
E = E / max(E(:));

IND = find(I);
[ROWS, COLS] = ind2sub([m, n], IND);

RESULTS = [ROWS, COLS, E(IND), U(IND), V(IND)];
% Sort before saving to file, C++ reports the links for the sorted points
P = round(sortrows(P, [1, 2])); 

dlmwrite('data.txt', [size(RESULTS, 1), size(P, 1), m, n, MAX_FILLS, GAMMA], 'delimiter', ' ', 'precision', 10)
dlmwrite('data.txt', RESULTS, '-append', 'delimiter', ' ', 'precision', 5)
dlmwrite('data.txt', P, '-append', 'delimiter', ' ', 'precision', 5)


% fid = fopen('neuron.bin', 'w');
% fwrite(fid, [size(RESULTS, 1), size(P, 1), m, n], 'int16');
% RESULTS(:, 3:end) = RESULTS(:, 3:end) * 10000;
% fwrite(fid, int16(RESULTS), 'int16');
% fwrite(fid, int16(P), 'int16');

%% Dijkstra Shortest Path

N = size(P, 1);

A       = zeros(N);
PATHS   = cell(N);

fprintf(1, 'Runnig the C++ code\n');
!./a.out < data.txt > result.txt
fprintf(1, 'Reading the results\n');

fid = fopen('result.txt', 'r');
while ~feof(fid)
    str = fgetl(fid);
    C = textscan(str, '%u32');
    path = reshape(C{1}, 2, []);
    path = path';
    
    A(path(1, 1) + 1, path(1, 2) + 1) = 1;
    A(path(1, 2) + 1, path(1, 1) + 1) = 1;
    PATHS{path(1, 1) + 1, path(1, 2) + 1} = path(2:end, :);
    PATHS{path(1, 2) + 1, path(1, 1) + 1} = path(2:end, :);
    
end
fclose(fid);


% Q = parallel.pool.DataQueue;
% step = 1;
% h = waitbar(0, 'Tracing connections...');
% afterEach(Q, @nUpdateWaitBar);
% 
% 
% p = gcp();
% 
% stepSize = idivide(uint32(N), p.NumWorkers, 'ceil');
% 
% for idx = 1:p.NumWorkers
%   F(idx) = parfeval(@DijkstraNeuron, 3, I, E, U, V, P, (idx - 1) * stepSize + 1:min(N, idx * stepSize), Q, MAX_FILLS);
% end
% 
% for idx = 1:p.NumWorkers
%   
%     % fetchNext blocks until next results are available.
%   [~, AF, AEF, PATHSF] = fetchNext(F);
%   
%   IDX           = ~isinf(AF);
%   A(IDX)        = AF(IDX);
%   AE(IDX)       = AEF(IDX);
%   PATHS(IDX)    = PATHSF(IDX);
% end
% 
% close(h);

%% Show results
M = zeros(size(I));

for i = 1:size(PATHS, 1)
    for j = 1:size(PATHS, 2)
        
        path = PATHS{i, j};
        
        for k = 1:size(path, 1)
            M(path(k, 1), path(k, 2)) = 1.0;
        end
    end
end

figure;
imshow(cat(3, M, IMG, I))

end

