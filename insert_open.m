function new_row = insert_open(xval,yval,parent_index,hn,gn,fn,...
                                target_index,theta)
%Function to Populate the OPEN LIST
%OPEN LIST FORMAT
%--------------------------------------------------------------------------
%IS ON LIST 1/0 |X val |Y val |Parent index |h(n) |g(n)|f(n)| target_index
%former theta |
%-------------------------------------------------------------------------
%
%   Copyright 2009-2010 The MathWorks, Inc.
new_row = [1, xval, yval, parent_index, hn, gn, fn, target_index, theta];

end