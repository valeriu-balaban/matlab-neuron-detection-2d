image = imread('preprocessed-t0.tif');
image = image(4001:5000, 4001:5000);

%% parameters

RADIUS_LIST = 25:60;
EDGE_WIDTH = 3;
THETA_THRESHOLD = 0.6;

%% find neurons

[position, radius, BW] = FindNeurons(image, RADIUS_LIST, EDGE_WIDTH, THETA_THRESHOLD);

%% assign possibility and direction for synapse
% Parameters

ALPHA = -1/3;
SIGMA = 1;
GAMMA = 0.8;

% Transformation matrix for the eigenvalues
T = [1, ALPHA; ALPHA, 1];

% same function as img_as_float
image = double(image);
image = image ./ max(image(:));

% Hessian computation
gaussian_img = imgaussfilt(image, SIGMA, 'padding', 'replicate');

% second-order derivatives
[gx, gy] = gradient(gaussian_img);
[gxx, gxy] = gradient(gx);
[gyx, gyy] = gradient(gy);

% Compute line directions
[m, n] = size(image);
E_min  = zeros(m, n);
E_max  = zeros(m, n);
U      = zeros(m, n);
V      = zeros(m, n);

for i = 1:m
   for j = 1:n 
        Hessian = [gxx(i,j) gxy(i,j); gxy(i,j) gyy(i,j)];
        [eigen_vector, eigen_value] = eig(Hessian, 'vector');
        
        % Adjust eigenvalues
        eigen_value = abs(T * eigen_value);
        
        E_max(i, j) = max(eigen_value);
        [E_min(i, j), idx] = min(eigen_value);
        
        % eigenvector
        U(i, j) = eigen_vector(1, idx);
        V(i, j) = eigen_vector(2, idx);
   end
end
theta = atan2d(V, U);

%%
BW_thin = bwmorph(BW, 'thin', inf);
% figure; imshow(BW_thin);
% hold on;
% quiver(U .* BW_thin, V .* BW_thin);
%% Astar algorithm for synapse detection
% draw neurons
draw_circles(position, radius, image);

% label
for i = 1:length(radius)
    text(position(i,1),position(i,2),int2str(i),'FontSize',10,'Color','red');
end
% find
theta_thre = 15;
fill_gap = 20;

%init
num = size(position, 1);
connected = zeros(num, num);
synapse = cell(num, num);

parfor i = 1:num
    [connected(i, :), synapse(i, :)] = Astar_Synapse(BW_thin, i, radius, ...
                                    position, theta, fill_gap, theta_thre); 
end

% make sure connectivity is symmetrical and least length
for i = 1:num-1
    for j = i+1:num
        if connected(i, j) && connected(j, i) && connected(i, j) > connected(j, i)...
                || connected(j, i) && ~connected(i, j)
            connected(i, j) = connected(j, i);
            synapse{i,j} = flip(synapse{j, i}, 1);
        end    
    end
end

% compute width of synapses
breadth = synapse_breadth(synapse, BW);

% draw synapses
draw_synapse(connected, synapse, breadth);