function [output_data, time_pt] = moving_time_window_v2(data, window_slide, size_window, func)
% Input data must be time x channel
% output_data is time window averaged data

% edited by Youngjai at July 28th, 2023.

% if you input the desirable function, you can get the functioned
% moving and globally averaged data as output.

% INPUT
% data         : time by channel, what to do moving averaged
% window_slide : the interval between two points, unit : # of points
% size_window  : the window size, unit : # of points
% func         : the desired function when to slice the data

% OUTPUT
% output_data  : rearranged time by channel
% time_pt      : reorganized time information

    size_data = size(data,1);
    size_ch = size(data,2);
    
    mid_time_pt = (size_window-1)/2+1;
    
    reindexes = 0:window_slide:size_data-size_window;
    time_pt = zeros( length(reindexes), 1);
    output_data = zeros( length(reindexes), size_ch);
    
    index = 1;
    if size_window ~=1
         for i=reindexes 
             time_pt(index) = mid_time_pt;
             output_data(index,:) = mean( func(data(i+1:size_window+i,:)), 1);
             mid_time_pt = mid_time_pt + window_slide;
             index = index + 1;
         end
    else
       output_data = func(data);
       time_pt = 1:size_data;
    end

end