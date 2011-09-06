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
% Starburst Algorithm - Version 1.0.0
% Part of the openEyes ToolKit -- http://hcvl.hci.iastate.edu/openEyes
% Release Date:
% Authors : Dongheng Li <donghengli@gmail.com>
%           Derrick Parkhurst <derrick.parkhurst@hcvl.hci.iastate.edu>
% Copyright (c) 2005
% All Rights Reserved.

function calculate_limbus_and_gaze()

% This function calculates the elliptical limbus parameters in the eye
% image and the preditive locacation of gaze in the scene image.
%
% The user need to specify the directory which contains the /Eye and /Scene image directory and the calibration
% file (calibration.mat). The user also need to eliminate bad calibration point correspondences  by click at the blue
% cross indicating the location of the corneal reflection). And then the program will calculate the gaze location
% in the scene image based.

pname = uigetdir(pwd,'Select Dir of data');
scene_file = sprintf('%s/Scene/Scene_', pname);
eye_file = sprintf('%s/Eye/Eye_', pname);
calibration_data_name = sprintf('%s/calibration.mat', pname);
ellipse_result_file = sprintf('%s/ellipse_result.mat', pname);
dat_file = sprintf('%s/ellipse_result.dat', pname);
eval(sprintf('!mkdir %s/Result_Small_Eye', pname));
small_eye_result_file = sprintf('%s/Result_Small_Eye/result_', pname);
libmus_edge_thresh = 15;
dis_step = 10;
start_frame = 205;

% Automatically get the frame number range
first_frame = get_first_or_last_frame_num(sprintf('%s/Eye/', pname), 'Eye_', 5, 'first');
last_frame = get_first_or_last_frame_num(sprintf('%s/Eye/', pname), 'Eye_', 5, 'last');
synch_eye_minus_scene = 3;  %synchronize eye and scene image
last_frame = last_frame - abs(synch_eye_minus_scene);

% Eliminate the bad calibration points
load(sprintf('%s', calibration_data_name));
bad_cal_indices = [];
Is = read_image(scene_file, cal_frame(end));
Ie = read_gray_image(eye_file, cal_frame(end)+synch_eye_minus_scene);
figure, 
subplot(1,2,1); imshow(uint8(Is)); hold on;
title({'Scene image:'; 'green cross: preditive gaze location'; 'Eye image:'; 'red cross: pupil center';
       'yellow star(*): bad calibration correspondence (if any)'});
plot(cal_scene(:,1), cal_scene(:,2), 'g+');
while 1,
    subplot(1,2,2); imshow(uint8(Ie)); hold on;
    title(sprintf(' Instruction:\n left click on red blue cross to eliminate bad \n calibration point correspondence\n when you finish, click other button'));
    plot(cal_ellipse(:,3), cal_ellipse(:,4), 'r+');
    if ~isempty(bad_cal_indices)
        plot(cal_ellipse(bad_cal_indices,1), cal_ellipse(bad_cal_indices,2), 'y*');
    end
    [tx,ty,but] = ginput(1);
    if but == 1,
       dis = sqrt((cal_ellipse(:,1)-tx).^2 + (cal_ellipse(:,2)-ty).^2);
       min_dis_index = find(dis==min(dis), 1, 'first');
       bad_cal_indices = [bad_cal_indices min_dis_index];
    else
        break;
    end
end
bad_cal_indices = [bad_cal_indices 0]; % add a zero in order to eliminate the case of empty set
    
% while 1,% Starburst Algorithm
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

%     n = input('input the index of bad calibration points obtained by looking at the calibration.jpg (should be 1-9; 0-end input)\n');
%     if n >= 1 && n <= 9
%        bad_cal_indices = [bad_cal_indices n];
%     elseif n == 0
%         break;
%     else
%         fprintf('Error input! should be 1-9 or 0 to finish\n');
%     end
% end

% % Calculate the homography mapping matrix
% % Use the different vector between pupil center and corneal reflection
% [neye_x, neye_y, T1] = normalize_point_coordinates(cal_ellipse(:,3)-cal_cr(:,1), cal_ellipse(:,4)-cal_cr(:,2));
% [ncal_x, ncal_y, T2] = normalize_point_coordinates(cal_scene(:,1), cal_scene(:,2));
% A = zeros(2*length(ncal_x), 9);
% for i=1:length(ncal_x),
%     if i ~= bad_cal_indices
%         A(i*2-1,:) = [0 0 0 -neye_x(i) -neye_y(i) -1 ncal_y(i)*neye_x(i) ncal_y(i)*neye_y(i) ncal_y(i)];
%         A(i*2,:) = [neye_x(i) neye_y(i) 1 0 0 0 -ncal_x(i)*neye_x(i) -ncal_x(i)*neye_y(i) -ncal_x(i)];
%     end
% end
% [ua, sa, va] = svd(A);
% c = va(:,end);
% H_cr=reshape(c,[3,3])';
% H_cr=inv(T2)*H_cr*T1

% Calculate the second order polynomial parameters for the mapping
% Use the different vector between pupil center and corneal reflection
eye_x = cal_ellipse(:,3);
eye_y = cal_ellipse(:,4);
cal_x = cal_scene(:,1);
cal_y = cal_scene(:,2);
A = zeros(length(cal_x), 6);
for i=1:length(cal_x),
    if i ~= bad_cal_indices
        A(i,:) = [eye_y(i)^2 eye_x(i)^2 eye_y(i)*eye_x(i) eye_y(i) eye_x(i) 1];
    end
end
[ua, da, va] = svd(A);
b1 = ua'*cal_x;
b1 = b1(1:6);
par_x = va*(b1./diag(da));
b2 = ua'*cal_y;
b2 = b2(1:6);
par_y = va*(b2./diag(da));

frame_index = start_frame;
Ie = read_gray_image(eye_file, frame_index);
fig_handle = figure, imshow(uint8(Ie));
title(sprintf('Please click near the limbus center'));
[cx, cy] = ginput(1);
close(fig_handle);

[height width] = size(Ie);
small_eye_ratio = 0.25;
small_h = uint16(height*small_eye_ratio);
small_w = uint16(width*small_eye_ratio);

scene = zeros(last_frame, 2);
ellipse = zeros(last_frame, 5);
small_eye_ratio = 0.25;
cross_len = 5;
red = [255 0 0];
green = [0 255 0];
blue = [0 0 255];
area = pi*cal_ellipse(:,1).*cal_ellipse(:,2);
area_avg = mean(area);
area_std = std(area);
ellipse_count = 9;
avg_num = 5;
scene_array = zeros(avg_num,2);
tic
for frame_index=start_frame:last_frame
    fprintf(1, '%d-', frame_index);
    if (mod(frame_index,30) == 0)
        fprintf(1, '\n');
    end

    Ie = read_gray_image(eye_file, frame_index+synch_eye_minus_scene);
    Ie_small = imresize(Ie, small_eye_ratio);
    Ie_small = repmat(Ie_small, [1 1 3]);
    Is = read_image(scene_file, frame_index);
    Ie = gaussian_smooth_image(Ie, 1.5);
    [ellipse(frame_index,:)] = detect_limbus_two_features(Ie, libmus_edge_thresh, cx, cy, -pi/5, pi/5, dis_step, area_avg, 1.5*area_std);
    
    if ~(ellipse(frame_index,1) <= 0 || ellipse(frame_index, 2) <= 0)
        cx = ellipse(frame_index, 3);
        cy = ellipse(frame_index, 4);
        
        % Calbulate the scene position using the pupil center
        coef_vector = [cy^2 cx^2 cy*cx cy cx 1];
        scene(frame_index,:) = [coef_vector*par_x coef_vector*par_y];
        %scene_pos = H_cr*[cx-cr(frame_index,1) cy-cr(frame_index,2) 1]';
        %scene(frame_index,:) = [scene_pos(1)/scene_pos(3)
        %scene_pos(2)/scene_pos(3)];
        
        if ellipse(frame_index,3) ~= 0 & ellipse(frame_index,4) ~= 0
            Ie_small = plot_ellipse_in_image(Ie_small, ellipse(frame_index,:)*small_eye_ratio, 'g', 2);
        end
        count = 0;
        for i=0:avg_num-1
            if scene(frame_index-i,1) > 0 & scene(frame_index-i,1) <= width & scene(frame_index-i,2) > 0 & scene(frame_index-i,2) <= height
                count = count+1;
                scene_array(count,:) = scene(frame_index-i,:);
            end
        end
        if count == 1
            scene_pos = round(scene_array(1,:));
        else
            scene_pos = round(mean(scene_array(1:count,:)));
        end
        if scene_pos(1) >= cross_len+2 & scene_pos(1) <= width-cross_len-1 & ...
                scene_pos(2) >= cross_len+2 & scene_pos(2) <= height-cross_len-1
            for bit = 1:3
                Is(scene_pos(2), scene_pos(1)-cross_len:scene_pos(1)+cross_len,bit) = green(bit);
                Is(scene_pos(2)-cross_len:scene_pos(2)+cross_len, scene_pos(1),bit) = green(bit);
                Is(scene_pos(2)-1,scene_pos(1)-cross_len:scene_pos(1)+cross_len,bit) = green(bit);
                Is(scene_pos(2)-cross_len:scene_pos(2)+cross_len, scene_pos(1)-1,bit) = green(bit);
                Is(scene_pos(2)+1,scene_pos(1)-cross_len:scene_pos(1)+cross_len,bit) = green(bit);
                Is(scene_pos(2)-cross_len:scene_pos(2)+cross_len, scene_pos(1)+1,bit) = green(bit);
            end
        end
    end
     
    Is(1:small_h, 1:small_w, :) = Ie_small(1:small_h, 1:small_w, :);
    imwrite(uint8(Is), sprintf('%s%05d.jpg', small_eye_result_file, frame_index));
end
toc
save(sprintf('%s', ellipse_result_file), 'ellipse', 'scene');
save_eye_tracking_data(dat_file, ellipse, scene);
for i=start_frame:last_frame, 
    eval(sprintf('!mv %s%05d.jpg %s%05d.jpg\n',small_eye_result_file,i,small_eye_result_file,i-start_frame)); 
end
eval(sprintf('!ffmpeg -i %s%%05d.jpg -b 1200 %s/small_eye_result.mpg', small_eye_result_file, pname));

function [I] = read_gray_image(file, index);
I = double(rgb2gray(imread(sprintf('%s%05d.jpg', file, index))));

function [I] = read_image(file, index);
I = double(imread(sprintf('%s%05d.jpg', file, index)));
