function [dist, index] = mydistance(current, Target)
%This function calculates the distance between current point and all target
%neurons and find the target with mimimum distance.

num = size(Target, 1);

[dist, index] = min(sqrt(sum((repmat(current, num, 1) - Target) .^ 2, 2)));

