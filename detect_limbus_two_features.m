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

function limbus_ellipse = detect_limbus_two_features(I, edge_thresh, cx, cy, start_angle, end_angle, dis_step, area_avg, area_delta)

% This function detects pupil and corneal reflection in the eye image
%
% Input:
% I = input image
% edge_thresh = threshold for pupil edge detection
% cx = x coordinate of starting point
% cy = y coordinate of starting point
% start_angle = the starting angle to detect features
% end_angle = the ending angle to detect features
% dis_step = distance step to detect feature
% area_avg = average area of the ellipse
% area_delta = range centered at the average area of the valid ellipse
%
% Output:
% limbus_ellipse = 5-vector of the ellipse parameters of pupil
%   [a b cx cy theta]
%   a - the ellipse axis of x direction
%   b - the ellipse axis of y direction
%   cx - the x coordinate of ellipse center
%   cy - the y coordinate of ellipse center
%   theta - the orientation of ellipse

angle_delta = 1*pi/180;         % discretization step size (radians)
min_feature_candidates = 10;    % minimum number of pupil feature candidates
max_ransac_iterations = 10000;  % maximum number of ransac iterations

epx = []; epy = [];
angle_step = pi/360;
[epx, epy] = limbus_feature_detection(I, cx, cy, edge_thresh, angle_step, start_angle, end_angle, dis_step);
[epx1, epy1] = limbus_feature_detection(I, cx, cy, edge_thresh, angle_step, pi-end_angle, pi-start_angle, dis_step);
epx = [epx epx1];
epy = [epy epy1];

if isempty(epx) || isempty(epy)
    limbus_ellipse = [0 0 0 0 0]';
    return;
end

[epx, epy] = distance_filter(epx, epy, 1.5, cx, cy);
[epx, epy] = distance_filter(epx, epy, 1.5, mean(epx), mean(epy));

if exist('area_avg') & exist('area_delta')
    [ellipse, inliers] = fit_ellipse_ransac_limbus(epx, epy, max_ransac_iterations, area_avg, area_delta);
else
    [ellipse, inliers] = fit_ellipse_ransac_limbus(epx, epy, max_ransac_iterations);
end
     
if ellipse(3) < 1 || ellipse(3) > size(I,2) || ellipse(4) < 1 || ellipse(4) > size(I,1)
    fprintf(1, 'Error! The ellipse center lies out of the image\n');
    limbus_ellipse = [0 0 0 0 0]';
else
    limbus_ellipse = ellipse;
end
