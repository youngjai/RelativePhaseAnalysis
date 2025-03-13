%script for quick plotting

%parameters
MAT = dpli_out; %assign the countourr matrix from topoplot_general
load('humanTopBorder.mat'); %load the appropriate border coordinates

resolution = size(MAT,1);
%scale
X = borderCoords(2,:);
Y = borderCoords(1,:);
Xlin = linspace(min(X), max(X), resolution); 
Ylin = linspace(min(Y), max(Y), resolution); 
[xx,yy] = meshgrid(Xlin, Ylin); 

colormap('jet');
contourf(xx,yy,MAT)

hold on; 
plot3(borderCoords(2,:),borderCoords(1,:), 100*ones(size(borderCoords,2)),'color',[0 0 0],'LineWidth',2);

axis off; 
axis equal; 