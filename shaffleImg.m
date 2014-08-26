function img = shaffleImg(img, mask, windowSize)

[row col] = find(mask == 1);

samplingImg = [repmat(img(:,1), [1,windowSize]), img, repmat(img(:,end), [1,windowSize])];
samplingImg = [repmat(samplingImg(1,:), [windowSize,1]); samplingImg; repmat(samplingImg(end,:), [windowSize,1])];  

for i=1:1:length(row)
    img(row(i), col(i)) = samplingImg(windowSize + row(i) + randi([-windowSize windowSize]), windowSize + col(i) + randi([-windowSize windowSize]));
end

end