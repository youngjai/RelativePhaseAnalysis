%open R05data_rel_rs.mat

time_moving=1000;
time_window=2000;
fdat_data = data2_30';
[Fdatm_phase, Fdatm_rel_p, Fdatm_rel_p_ave, Fdatm_phase_ave] = rel_phase_cal_mov( fdat_data, time_moving, time_window);

rel_p =  double(Fdatm_phase_ave);

k=1;

vidfile = VideoWriter([ ' brain_data2_30_v' num2str(k) ] ,'MPEG-4');
vidfile.FrameRate = round(1000/time_moving);
open(vidfile);    

time_dt = 1/1000*time_moving
time_now = 0;

%for i=1:size(rel_p,1)
 dt=1;
 tic
for i=1:dt:size(rel_p,1)
    
figure(1);
topoplot_general_test(rel_p(i,:), chanloc_cord,'smooth',25);
title(  [ '\rm' num2str(time_now) 'S'; ] ,'fontsize',16 )
drawnow
F(i) = getframe(gcf); 
writeVideo(vidfile,F(i));
time_now = time_now + dt*time_dt;
end
 close(vidfile)
toc