%flip coord on the x plane
minX = min(borderCoords(1,:));
borderCoords(1,:) = -borderCoords(1,:);
difference = minX - min(borderCoords(1,:));
borderCoords(1,:) = borderCoords(1,:) + difference;