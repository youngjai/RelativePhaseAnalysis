function inner_border = generate_lining(borderCoords, percent)
% inner_border: Generates a set of coordinates that will line the inner
%               side of the given shape.  For example, given the
%               coordinates for a simple square, inner_border will generate
%               coordinates of a smaller square that will fit right inside
%               the original square.
%
%   borderCoords must be a set of ordered coordinates in which plotting
%   them will result in a closed polygon.  
inner_border = borderCoords;

%find center
maxX = max(borderCoords(1,:));
minX = min(borderCoords(1,:));
maxY = max(borderCoords(2,:));
minY = min(borderCoords(2,:));

centerx = minX + (maxX - minX) / 2;
centery = minY + (maxY - minY) / 2;

%shrink border coords towards the center
for i = 1:size(borderCoords,2)
    if(borderCoords(1,i) > centerx)
        inner_border(1,i) = inner_border(1,i) - (inner_border(1,i) - centerx) * percent;
    else
        inner_border(1,i) = inner_border(1,i) + (centerx - inner_border(1,i)) * percent;
    end
    
    if(borderCoords(2,i) > centery)
        inner_border(2,i) = inner_border(2,i) - (inner_border(2,i) - centery) * percent;
    else
        inner_border(2,i) = inner_border(2,i) + (centery - inner_border(2,i)) * percent;
    end
end

            
        