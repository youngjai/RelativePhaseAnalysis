function coords = generate_circle(x0, y0, r)
%generate_circle(x0, y0, r)
%   Description:    returns a list of coordinates for a circle centered on x0
%                   and y0 with radius r
angle=-pi:0.1:pi;
angl=angle(randi(numel(angle),100,1));
coords(2,:)=r*cos(angl)+x0;
coords(1,:)=r*sin(angl)+y0;