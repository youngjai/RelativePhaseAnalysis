function topoplot_figure(topo, borderCoords, xx, yy, Coord, varargin)
%topoplot_figure(topo, varagin)
%   Description:    Topological plot modified for general use with various
%                   shapes
% 
%   REQUIRED ARGUMENTS
%   topo:           Topoplot results in 200x200 format
%
%   xx, yy:         xx and yy coordinates from topoplot_general
%   
%   borderCoords:   Coordinates of the border
%
%   Coord:          Coordinates of the channels
%
%   ADDITIONAL PARAMETERS
%
%   border:         Toggles whether to plot the border or not.  1 = on, 
%                   0 = off.  Default is 1
%
%   scatter:        Toggles whether to plot scatter of coordinates 1 = on,
%                   0 = off.  Default is 0
%
%   color:          Set's the color scale for the topoplot.  0 = auto (relative),
%                   1 = absolute small range, 2 = absolute large range (true to regional plot)
%
%   line:           Toggles the border lines between fills of the contour
%   
p = inputParser;
p.addRequired('topo', @(x) length(x) > 0);
p.addRequired('xx', @(x) length(x) > 0);
p.addRequired('yy', @(x) length(x) > 0);
p.addRequired('Coord', @(x) length(x) > 0);
p.addRequired('borderCoords', @(x) length(x) > 0);

defaultborder = 1;
p.addParameter('border', defaultborder, @isnumeric);
defaultscatter = 0;
p.addParameter('scatter', defaultscatter, @isnumeric);
defaultColor = 0;
p.addParameter('color', defaultColor, @isnumeric);
defaultLine = 0;
p.addParameter('line', defaultLine, @isnumeric);

p.parse(topo, borderCoords, xx, yy, Coord, varargin{:});

border = p.Results.border;
scatterFlag = p.Results.scatter;
colorFlag = p.Results.color;
lineFlag = p.Results.line;

colormap(jet(128));
[hC hC] = contourf(xx,yy,topo); 
if (lineFlag == 0)
    set(hC,'LineStyle', 'none');
end
%color scale
if(colorFlag == 1)
    minValue = min(abs([max(max(topo)), min(min(topo))]));
    caxis([-minValue, minValue]);
elseif(colorFlag == 2)
    maxValue = max(max(abs(topo)));
    caxis([-maxValue,maxValue]);
end

%if custom color limits are supplied
% if clim(1) ~= -1 & clim(2) ~= 1
%     caxis(clim);
% end

 %colorbar('SouthOutside'); 
hold on; 
if border == 1
    plot3(borderCoords(2,:),borderCoords(1,:), 100*ones(size(borderCoords,2)),'color',[0 0 0],'LineWidth',2);
end

if scatterFlag == 1
    scatter(Coord(:,1), Coord(:,2),25, [0 0 0], 'filled');
end

axis off; 
axis equal; 
hold off; 