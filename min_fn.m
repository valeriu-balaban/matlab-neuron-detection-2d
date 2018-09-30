function [min_point_index, target_index] = min_fn(OPEN)
%Function to return the Node with minimum fn
% This function takes the list OPEN as its input and returns the index of the
% node that has the least cost
%

 OPEN_valid = OPEN(:, 1) == 1;
 
 %Send the index of the smallest node
 if sum(OPEN_valid)
     [~, min_index] = min(OPEN(OPEN_valid, 7));
     valid_index = find(OPEN_valid);
     min_point_index = valid_index(min_index);
     target_index = OPEN(valid_index, 8);
 else
     min_point_index = -1;
     target_index = -1;
 end