function dataOut = topoplot_general_test(temp_Z, Coord, varargin) 
%topoplot_general(temp_Z, Coord, borderCoords, varagin)
%   Description:    Topological plot modified for general use with various
%                   shapes
% 
%   REQUIRED ARGUMENTS
%   temp_Z:         values coresponding to each coordinate.  Must be same
%                   length as coordinates
%
%   Coord:          XY coordinates.  Supply in two columns (Matrix of size Nx2)
%                   must be the same length as temp_Z
%
%   borderCoords:   Coordinates for the outline of the animal.  Coodinates
%                   must be ordered such that plotting the coordinates will 
%                   form a closed shape.  Additionally, complex border shapes
%                   may not be filled correctly.  Several default coordinates have 
%                   been supplied.  Supply as two rows (Matrix of size 2xN)
%   
%   OPTIONAL ARGUMENTS
%   clim:           Range of values.  Default is [-1,1];
%
%   color_arg:      Assorted color schemes.  Default is Jet 
%                   1 = Jet, 2 = grey, 3 = autumn   
%
%   label_arg:      Names of the coodinates. Default is simply the integers 1:N
%
%   ADDITIONAL PARAMETERS
%   smooth:   Radius supplied to reduce resolution function.  
%                   The larger the number, the smoother the contours will
%                   be. Default is 1 (no smoothing)
%
%   expand:         supply 0 or 1 or 2.  1 expands to the border. 2 enables 
%                   the topoplot to be expanded to a bigger circle. Defaults to 0
%
%   shiftpreset:    shifts that have set to default borders.  Current options are
%                   1 = macaque, 2 = mouse, 3 = human top, 4 = human side,
%                   5 = human theory top, 6 = human side squeezed, 7 = human top squeezed, 8 = mouse theory
%                   9 = monkey theory, 10 = 96ch ketamine, 11 = human side simulation 12 = human theory side reconfigured
%                   13 = human experimental side reconfigured, 14 = mouse
%                   148 ch
%                   
%                   NOTE: will only work properly with presupplied border
%                   coordinates.
%
%   edge:           Edge correction trims the areas of the topoplot that
%                   bleed out of the border. 1 = on, 0 = off.  Default is 1 
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
%   extrapolate:    Toggles method of determining value of border
%                   coordinates. 1 = extrapolates to nearest coordinate. 0
%                   = average value of all non-border coordinates. Default
%                   is 0
%   
%   res:            Set's resolution of interpolated area.  Default is
%                   200x200
%
%   line:           Toggles the border lines between fills of the contour
%   
%
% Feb 2, 2016 by Joseph Lee

%% Parse Inputs
p = inputParser;
p.addRequired('temp_Z', @(x) length(x) > 0);
p.addRequired('Coord', @(x) length(x) > 0);
% p.addRequired('borderCoords', @(x) length(x) > 0);

defaultclim = [-1 1];
p.addOptional('clim', defaultclim, @isnumeric);

defaultcolor = 1;
p.addOptional('color_arg', defaultcolor, @isnumeric);
%default label.  Will adjust to right size later
defaultlabel = 1;
p.addOptional('label_arg', defaultlabel, @isnumeric);

%extra parameters
defaultradius = 1;
p.addParameter('smooth', defaultradius, @isnumeric);
defaultexpand = 1;
p.addParameter('expand', defaultexpand, @isnumeric);
defaultshift = [0 0 0 0];
p.addParameter('shift', defaultshift, @isnumeric);
defaultshiftpreset = 0;
p.addParameter('shiftpreset', defaultshiftpreset, @isnumeric);
defaultedge = 1;
p.addParameter('edge', defaultedge, @isnumeric);
defaultborder = 1;
p.addParameter('border', defaultborder, @isnumeric);
defaultscatter = 0;
p.addParameter('scatter', defaultscatter, @isnumeric);
defaultExtrapolate = 0;
p.addParameter('extrapolate', defaultExtrapolate, @isnumeric);
defaultColor = 0;
p.addParameter('color', defaultColor, @isnumeric);
defaultRes = 200;
p.addParameter('res', defaultRes, @isnumeric);
defaultLine = 0;
p.addParameter('line', defaultLine, @isnumeric);

%parse and fix defaults
p.parse(temp_Z, Coord, varargin{:});

clim = p.Results.clim;
color_arg = p.Results.color_arg;
label_arg = p.Results.label_arg;
resolutionRadius = p.Results.smooth;
expand = p.Results.expand;
shift = p.Results.shift;
shiftpreset = p.Results.shiftpreset;
edge = p.Results.edge;
border = p.Results.border;
scatterFlag = p.Results.scatter;
extrapolateFlag = p.Results.extrapolate;
colorFlag = p.Results.color;

resolution = p.Results.res;
lineFlag = p.Results.line;

if (isempty(label_arg))
    label = num2str((1:length(Coord))');
end

%parse shifting preset;
if shiftpreset == 0
    %no shift
   % shift = [0 0 0 0];
   load('borderCoords.mat')
elseif shiftpreset == 1
    %monkey
    shift = [-50 50 -60 100];
    load('monkeyBorder.mat');
elseif shiftpreset == 2
    %mouse
    shift = [-35 35 -75 100];
    load('mouseBorderInner.mat')
elseif shiftpreset == 3
    %human top
    shift = [25 -25 25 -25];
    load('humanTopBorder.mat')
elseif shiftpreset == 4
    %human side
    shift = [25 -25 25 -25];
    load('humanSideBorder.mat')
elseif shiftpreset == 5
    %human theory top
    shift = [-30 30 -50 50];
    load('humanTopBorder.mat')
elseif shiftpreset == 6
    %human side squeezed to fit into brain
    shift = [-50 50 -25 125];
    % shift = [-50 50 -25 50];
    load('humanSideBorder.mat')
elseif shiftpreset == 7
    %human top squeezed
    shift = [-25 25 -25 25];
    load('humanTopBorder.mat')
elseif shiftpreset == 8
    %mouse theory
    shift = [-100 100 -120 120];
    load('mouseBorderInner.mat');
elseif shiftpreset == 9
    %monkey theory
    shift = [-80 40 -75 45];
    load('monkeyBorder.mat')
elseif shiftpreset == 10
    %96ch ketamine shift
    shift = [-25 25 -25 25];
    load('humanTopBorder.mat')
elseif shiftpreset == 11
    %human side simulation to match channel plot
    shift = [-80 60 -15 55];
    load('humanSideBorder.mat');
elseif shiftpreset == 12
    %human side border reconfigured
    shift = [-80 60 -25 75];
    load('humanSideBorderTall.mat');
elseif shiftpreset == 13
    %human experimental side reconfiguredd
    shift = [-50 50 -25 145];
    load('humanSideBorderTall.mat');
elseif shiftpreset == 14
    %mouse 148ch
    shift = [-25 25 -50 50];
    load('mouseBorderInner.mat')
end
%parse color
if color_arg == 0 
    black = [0,0,0;.5, .5, .5; 1,1,1]; 
    colormap(black);    
elseif color_arg == 1
    colormap(jet(128));
elseif color_arg == 2
    colormap(gray(128));
elseif color_arg == 3
    colormap(autumn(128));
end

if isempty(clim) 
    clim = [fix(min(min(temp_Z))/10)*10, fix(max(max(temp_Z))/10)*10]
end

%% load boundary and generate inner boundary
L = length(Coord); 
inner_bound = generate_lining(borderCoords, 0.05);
x2 = inner_bound(2,:);
y2 = inner_bound(1,:);
L2 = length(inner_bound);

%normalize the data
bord_maxX = max(borderCoords(2,:));
bord_minX = min(borderCoords(2,:));
bord_maxY = max(borderCoords(1,:));
bord_minY = min(borderCoords(1,:));

maxXnew = bord_maxX + shift(1);
minXnew = bord_minX + shift(2);
maxYnew = bord_maxY + shift(3);
minYnew = bord_minY + shift(4);

maxX = max(Coord(:,1));
minX = min(Coord(:,1));
maxY = max(Coord(:,2));
minY = min(Coord(:,2));

Coord(:,1) = ((maxXnew - minXnew) / (maxX - minX)) * (Coord(:,1) - maxX) + maxXnew;
Coord(:,2) = ((maxYnew - minYnew) / (maxY - minY)) * (Coord(:,2) - maxY) + maxYnew;

%% mesh the space
X = Coord(:,1); 
Y = Coord(:,2); 

%expand topoplot to border
if(expand == 1)
    X(L+1:L+length(borderCoords)) = borderCoords(2,:);
    Y(L+1:L+length(borderCoords)) = borderCoords(1,:);
    
    %dummy
    x2(L2+1:L2+length(borderCoords)) = borderCoords(2,:);
    y2(L2+1:L2+length(borderCoords)) = borderCoords(1,:);
elseif(expand == 2)
    %expand topoplot to larger circle
    radius = 0.6 * max([bord_maxX - bord_minX, bord_maxY - bord_minY]);

    centerx = bord_minX + ((bord_maxX - bord_minX) / 2);
    centery = bord_minY + ((bord_maxY - bord_minY) / 2);
    coordCircle = generate_circle(centerx, centery, radius);
    X(L+1:L+length(coordCircle)) = coordCircle(2,:);
    Y(L+1:L+length(coordCircle)) = coordCircle(1,:);
end

Xlin = linspace(min(X), max(X), resolution); 
Ylin = linspace(min(Y), max(Y), resolution); 
[xx,yy] = meshgrid(Xlin, Ylin); 
a = 1:L;   a=a';
Z = temp_Z; 

%dummy plot
Xlin2 = linspace(min(x2), max(x2), resolution); 
Ylin2 = linspace(min(y2), max(y2), resolution); 
[xx2,yy2] = meshgrid(Xlin2, Ylin2);
z2 = ones(L2,1);

%fill extra coordinates with dummy value
if (expand == 2)
    %for larger circle
    Z(L+1:L+length(coordCircle)) = min(Z) + (max(Z) - min(Z)) / 2;
elseif (expand == 1)
    %for border
%     Z(L+1:L+length(borderCoords)) = min(Z);


     % extrapolate values for the border
     if (extrapolateFlag == 1)
        for i = 1:length(borderCoords)
            Z(L+i) = extrapolate(temp_Z, Coord, borderCoords(:,i));
        end
     else
         Z(L+1:L+length(borderCoords)) = min(Z) + (max(Z) - min(Z)) / 2;
     end
     
     z2(L2+1:L2+length(borderCoords)) = 100;
end

% prctile(temp_Z,5); 
z = griddata(X,Y,Z, xx,yy,'cubic'); 


%compare dummy grid and actual grid to clean up edges
if(edge == 1)
    %dummy grid for triming edges
    z2 = griddata(x2',y2', z2, xx2, yy2, 'cubic');
    for i = 1:resolution
        for j = 1:resolution
            if(z2(i,j) >= 95)
                z(i,j) = NaN;
            end
        end
    end
end

% clean up stray elements
if(edge == 1)
    %dummy grid for triming edges
    for i = 2:resolution-1
        for j = 2:resolution-1
            if(~isnan(z(i,j)) && isnan(z(i+1,j)) && isnan(z(i-1,j)))
                z(i,j) = NaN;
            end
        end
    end
end

%smooth
% z_trans = reduce_resolution_circle_unweighted(z,resolutionRadius,0.2);
z_trans = reduce_resolution_circle_singlethread(z,resolutionRadius,0.2);
dataOut = z_trans;
%graph
[hC hC] = contourf(xx,yy,z_trans); 
if (lineFlag == 0)
    set(hC,'LineStyle', 'none');
end
%color scale
if(colorFlag == 1)
    minValue = min(abs([max(max(z_trans)), min(min(z_trans))]));
    caxis([-minValue, minValue]);
elseif(colorFlag == 2)
    maxValue = max(max(abs(z_trans)));
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
% scatter(X, Y);
% if color_arg == 0 
%     text(Coord(:,1) - 0.2, Coord(:,2) - 0.3, label,'Color','r'); 
% end
% text(-0.2,0.2,'{\bf\it BP}'); 
% text(-0.2,-3.8,'{\bf\it LP}'); 
axis off; 
axis equal; 
hold off; 



