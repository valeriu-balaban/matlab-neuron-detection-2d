function draw_synapse(connected, synapse, breadth)

num = length(connected);
hold on;
title('synapse');

for i = 1 : num-1
    for j = i+1 : num
        if connected(i, j)
            len = size(synapse{i, j}, 1);
            bread_path = breadth{i, j};
            path = synapse{i, j};
            for k = 1:len-1
                plot([path(k, 1), path(k+1, 1)],[path(k, 2), path(k+1, 2)],'LineWidth',2,'Color','green');
               % text(path(k, 1),path(k, 2),int2str(bread_path(k)),'FontSize',10,'Color','red');
            end
        end
    end
end
            