%
function value = extrapolate(data, coord, target)
%find distances to all points and pick the closest point
distances = sqrt((coord(:,1) - target(2)).^2 + (coord(:,2) - target(1)).^2);

%find smallest distance and find corresponding value
min = realmax;
for i = 1:length(coord)

    if(distances(i) < min)
        value = data(i);
        min = distances(i);
        index = i;
    end
end

end


    