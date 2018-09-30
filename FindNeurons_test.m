img = imread('neurons-step1-result.tif');
img = img(3001:4000, 4001:5000);

%% parameters

RADIUS_LIST = 25:60;
EDGE_WIDTH = 3;
THETA_THRESHOLD = 0.6;

%% find neurons

[position, radius, BW] = FindNeurons(img, RADIUS_LIST, EDGE_WIDTH, THETA_THRESHOLD);

%% assign possibility and direction for synapse
BW_thin = bwmorph(BW, 'thin', 20);

sigma = 8;
[U, V, theta] = neurite_vector(BW, sigma);

theta = theta .* BW_thin;

%% Astar algorithm for synapse detection
% draw neurons
draw_circles(position, radius, img);

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

% compute width of synapses
breadth = synapse_breadth(synapse, BW);

% draw synapses
draw_synapse(connected, synapse, breadth);

 %% Synapse detection
% % degree and distance threshold when connecting broken synapse
% theta_thre = 15;
% fill_gap = 0;
% [connected, synapse] = synapse_detection(BW, BW_thin, position,...
%                                                 radius, theta, theta_thre, fill_gap);
% 
% % compute width of synapses
% breadth = synapse_breadth(synapse, BW);