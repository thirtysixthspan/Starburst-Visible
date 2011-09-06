% Starburst Algorithm
%
% This source code is part of the starburst algorithm.
% Starburst algorithm is free; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% Starburst algorithm is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with cvEyeTracker; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%
% Starburst Algorithm for Visible Spectrum Eye Tracking - Version 1.0.0
% Part of the openEyes ToolKit -- http://hcvl.hci.iastate.edu/openEyes
% Release Date:
% Authors : Dongheng Li <donghengli@gmail.com>
%           Derrick Parkhurst <derrick.parkhurst@hcvl.hci.iastate.edu>
% Copyright (c) 2006
% All Rights Reserved.

function [nx, ny] = distance_filter(x, y, dev_factor, cx, cy)

% This function is to filter the points with a distance filter according to
% the deviation factor.
% Input:
% x = the x coordinates of points
% y = the y coordinates of points
% dev_factor = deviation factor; if set to 0.5, it means the valid points' 
% distance to the average center is (-0.5std, 0.5std); std: standard
% deviation
% Output:
% nx = the x coordinates of points after filtering
% ny = the y coordinates of points after filtering

%cx = mean(x);
%cy = mean(y);
dis = sqrt((x-cx).^2 + (y-cy).^2);
mean_dis = mean(dis);
std_dis = std(dis);
min_dis = mean_dis - dev_factor*std_dis;
max_dis = mean_dis + dev_factor*std_dis;
index = dis > min_dis & dis < max_dis;
nx = x(index);
ny = y(index);

