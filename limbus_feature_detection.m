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

function [epx, epy] = limbus_feature_detection(I, cx, cy, edge_thresh, angle_step, angle_start, angle_end, dis_step)
% Input
% I = input image
% [cx cy] = central start point of the feature detection process
% edge_threshold = threshold to detect limbus
% angle_step = angle step to detect features
% angle_start = starting angle to detect features
% angle_end = ending angle to detect features (e.g.: if angle_start=-PI/2 and angle_end=PI/2, then the program detect features 
%       in the angle [-PI/2, PI/2] with respect to the start_point
% dis_step = distance step to detect features
% Ouput
% epx = x coordinate of feature candidates [row vector] 
% epy = y coordinate of feature candidates [row vector] 

[height width] = size(I);
epx = [];
epy = [];
dir = [];
ep_num = 0;  % ep stands for edge point
for angle = angle_start : angle_step : angle_end
    dis_cos = dis_step*cos(angle);
    dis_sin = dis_step*sin(angle);
    pp_cos = cos(angle);
	pp_sin = sin(angle);
    p1 = [round(cx+dis_cos) round(cy+dis_sin)];
    if p1(2) > height | p1(2) < 1 | p1(1) > width | p1(1) < 1
        continue;
    end
    pixel_value1 = I(p1(2),p1(1));
    detected_num = 0;
    step = 2;
    while detected_num < 2, %(rand(1,1)>0.01)+1,
        p2 = [round(cx+step*dis_cos) round(cy+step*dis_sin)];
        if p2(2) > height | p2(2) < 1 | p2(1) > width | p2(1) < 1
            break;
        end        
        pixel_value2 = I(p2(2),p2(1));
        if (pixel_value2 - pixel_value1 > edge_thresh),
            ep_num = ep_num+1;
            p1 = [cx+(step-1)*dis_cos+pp_cos cy+(step-1)*dis_sin+pp_sin];
            cur_pixel = I(round(p1(2)), round(p1(1)));
            max_diff = cur_pixel - pixel_value1;
            tepx = p1(1)-pp_cos;    % edge point x coordinate
            tepy = p1(2)-pp_sin;    % edge point y coordinate
            for i = 2:dis_step
                pixel_value1 = cur_pixel;
                p1 = [p1(1)+pp_cos p1(2)+pp_sin];
                cur_pixel = I(round(p1(2)), round(p1(1)));
                if (max_diff < cur_pixel - pixel_value1)
					tepx = p1(1)-pp_cos;    % edge point x coordinate
                    tepy = p1(2)-pp_sin;    % edge point y coordinate
					max_diff = cur_pixel - pixel_value1;
                end
            end
            
            if isempty(find(round(tepx)==round(epx) &  round(tepy)==round(epy)))
                epx(ep_num) = tepx;
                epy(ep_num) = tepy;
            else
                ep_num = ep_num-1;
            end
            
            p1 = [round(p1(1)+2*dis_cos) round(p1(2)+2*dis_sin)];
            if p1(2) > height | p1(2) < 1 | p1(1) > width | p1(1) < 1
                break
            end
            pixel_value1 = I(p1(2),p1(1));
            step = step+2;
            detected_num = detected_num+1;
        else
            pixel_value1 = pixel_value2;
            step = step+1;
        end
    end
end
