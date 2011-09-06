function plot_limbus_and_avg_gaze()

clear all; close all;
pname = uigetdir(pwd,'Select Dir of data');
scene_file = sprintf('%s/Scene/Scene_', pname);
eye_file = sprintf('%s/Eye/Eye_', pname);
ellipse_result_file = sprintf('%s/ellipse_result.mat', pname);
load(ellipse_result_file);
eval(sprintf('!mkdir %s/Result_Small_Eye_Avg', pname));
small_eye_result_file = sprintf('%s/Result_Small_Eye_Avg/result_', pname);
libmus_edge_thresh = 15;

% Automatically get the frame number range
first_frame = get_first_or_last_frame_num(sprintf('%s/Eye/', pname), 'Eye_', 5, 'first');
last_frame = get_first_or_last_frame_num(sprintf('%s/Eye/', pname), 'Eye_', 5, 'last');
start_frame = 2620;
frame_index = start_frame;

cross_len = 5;
red = [255 0 0];
green = [0 255 0];
blue = [0 0 255];
Ie = read_gray_image(eye_file, frame_index);
[height width] = size(Ie);
small_eye_ratio = 0.25;
small_h = uint16(height*small_eye_ratio);
small_w = uint16(width*small_eye_ratio);
avg_num = 5;
scene_array = zeros(avg_num,2);
tic
for frame_index=start_frame:last_frame
    fprintf(1, '%d-', frame_index);
    if (mod(frame_index,30) == 0)
        fprintf(1, '\n');
    end

    Ie = read_gray_image(eye_file, frame_index);
    Ie_small = imresize(Ie, small_eye_ratio);
    Ie_small = repmat(Ie_small, [1 1 3]);
    Is = read_image(scene_file, frame_index);

    count = 0;
    if ~(ellipse(frame_index,1) <= 0 || ellipse(frame_index, 2) <= 0)
        for i=0:avg_num-1
            if scene(frame_index-i,1) > 0 & scene(frame_index-i,1) <= width & scene(frame_index-i,2) > 0 & scene(frame_index-i,2) <= height
                count = count+1;
                scene_array(count,:) = scene(frame_index-i,:);
            end
        end
        
        if ellipse(frame_index,3) ~= 0 & ellipse(frame_index,4) ~= 0
            Ie_small = plot_ellipse_in_image(Ie_small, ellipse(frame_index,:)*small_eye_ratio, 'g', 2);
        end
        scene_pos = round(mean(scene_array));
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
    imwrite(uint8(Is), sprintf('%s%05d.jpg', small_eye_result_file, frame_index-start_frame));
end
toc
eval(sprintf('!ffmpeg -i %s%%05d.jpg -b 400 %s/small_eye_result.mpg', small_eye_result_file, pname));


function [I] = read_gray_image(file, index);
I = double(rgb2gray(imread(sprintf('%s%05d.jpg', file, index))));

function [I] = read_image(file, index);
I = double(imread(sprintf('%s%05d.jpg', file, index)));
