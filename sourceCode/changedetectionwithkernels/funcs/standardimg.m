%  [matstd,mu,sd] = standardimg(matrix,t)
%  
%  IN:
%  matrix = image to normalize, 
%  t 	  = type of normalization (maxval  = pixel_value./max_pixel_value of the image)
%            			  (zscores = (x(band) - mean(band))/std(band))
%				  (zeroone = [0,1] \forall x_i)  
%				  (minusonneone = same as above, [-1,1])
% 				  (center  = x(band) - mean(
% OUT: matstd = normalized data
%      mu,sd     = mean and standard deviation (if used)
 
function [matstd,mu,sd] = standardimg(matrix,t)

[s1 s2 s3] = size(matrix);

if s3 > 1
    matrix = reshape(matrix,s1*s2,s3);
elseif s1 == 1
    matrix = matrix';
end

switch  t
    case 'maxval'   % data normalization ./ max val
        maxval = (max(max(matrix)));
        matstd = matrix./maxval;
    case 'zscores'   % data normalization
        [matstd,mu,sd] = zscore(matrix);
    case 'zeroone'   % scale data in [0,1]
        mx = max(matrix,[],1);
        mi = min(matrix,[],1);
        matstd = (matrix - repmat(mi,size(matrix,1),1)) ./ repmat(mx-mi,size(matrix,1),1);
    case 'minusoneone'    % scale data in [-1,1]
        for k = 1:s2;
            col = matrix(:,k);
            mx = max(col);
            mi = min(col);
            d = 2/(mx - mi);
            C = ( -mx - mi)/(mx - mi);
            redcol = col * d + C;
            matstd(:,k) = redcol;
        end
    case 'center'   % data normalization
        me = mean(matrix,1);
        matstd = matrix - repmat(me,size(matrix,1),1);
end

if s3 > 1
    matstd = reshape(matstd,s1,s2,s3);
elseif s1 == 1
    matstd = matstd';
end
