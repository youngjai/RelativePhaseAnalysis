% Input data must be time x channel
% output_data is time window averaged data

function [ output_data ,time_pt] = moving_time_window( data, window_slide, size_window)

time_rel_phase = size(data,1);
ch_rel_phase = size(data,2);

mid_time_pt = (size_window-1)/2;


time_pt = zeros( length(0:window_slide:time_rel_phase - size_window) , 1);
output_data = zeros( length(0:window_slide:time_rel_phase - size_window) , ch_rel_phase );

index = 1;

if size_window ~=1

     for i=0:window_slide:time_rel_phase - size_window
        
         time_pt(index) = mid_time_pt;
    
         output_data(index,:) = mean( data(i+1 : size_window + i , : ) ,1 );
        
         mid_time_pt = mid_time_pt + window_slide;
        
         index = index + 1;
    
     end

else
    output_data = data;
    time_pt = 0:1:time_rel_phase-1;
end