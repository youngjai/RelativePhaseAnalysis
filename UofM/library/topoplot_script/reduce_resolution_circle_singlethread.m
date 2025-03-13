function [z_trans ] = reduce_resolution_circle_singlethread(z,factor,epsilon)


size1z = size(z,1);
size2z = size(z,2);
z_trans = nan*ones(size1z,size2z);
%factor = 2;

z_bigger = [ nan*ones(size1z,factor) z nan*ones(size1z,factor) ] ;
z_bigger = [ nan*ones(factor,size1z+2*factor); z_bigger ; nan*ones(factor,size1z+2*factor) ] ; 

%generate disance matrix
weight_mat = zeros((factor * 2) + 1, (factor * 2) + 1);
for x=1:(factor * 2) + 1
    for y = 1:(factor * 2) + 1
        weight = sqrt((factor + 1 - x)^2 +(factor + 1 - y)^2);
        if(weight > factor)
            weight_mat(x,y) = nan;
        elseif(weight == 0)
            weight_mat(x,y) = 1;
        else
            weight_mat(x,y) = weight;
        end
    end
end

for i=1:size1z
    for j = 1:size2z
        if isnan (z(i,j)) == 0 ;
            z_temp = z_bigger( i+factor-factor : i +2*factor , j+factor-factor : j +2*factor ) ;
            z_nan = z_temp;
            z_nan(~isnan(z_temp)) = 1;
            %mean the diamond around center
            weightSum = 0;
            valueSum = 0;
            
            valueSum = nansum(nansum(z_temp ./ weight_mat));
            weightSum = nansum(nansum(1 ./ (weight_mat .* z_nan)));
%             for x=1:(factor * 2) + 1
%                 for y = 1:(factor * 2) + 1
%                     
%                     %check for nans
%                     if(isnan(z_temp(x,y)))
%                         continue;
%                     end
%                    
%                     %calculate distance
%                     weight = sqrt((factor + 1 - x)^2 +(factor + 1 - y)^2);
%                     
% %                     check weight for circle
%                     if(weight > factor)
%                         continue;
%                     end
%                               
%                     if(weight ~= 0)
%                         valueSum = valueSum + z_temp(x,y) / weight;
%                         weightSum = weightSum + (1 / weight);
%                     else
%                         valueSum = valueSum + z_temp(x,y);
%                         weightSum = weightSum + 1;
%                     end
%                     
%                end
%             end
            
            z_trans(i,j) = valueSum / weightSum;
        end
    end
end