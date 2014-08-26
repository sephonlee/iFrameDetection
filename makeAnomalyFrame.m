function synFrame = makeAnomalyFrame(strongModel, avgModel, numFrames, outVariables, needShaffle, randWindow)

if nargin == 5
    randWindow = floor(numFrames / 4);
elseif nargin == 4
    needShaffle = 0;
end

if needShaffle == 1
    fprintf('shaffle model...\n');
    strongModel.varImg = shaffleImg(strongModel.varImg, strongModel.maskActPixel, 5);
    fprintf('New model completed...\n');
end


startFrameId = outVariables.id + 1;
endFrameId = outVariables.id + numFrames;

isReverse = randi([1 3]);
if isReverse - 1 == 0 
    fprintf('Reverse Signal\n');
end


for ii = 1:1:numFrames
    
    if isReverse - 1 > 1;
        i = numFrames - ii;
    else
        i = ii;
    end
    
    j = i + randi([-randWindow randWindow]);
    if j > numFrames
        j = numFrames;
    elseif j < 0
        j = 0;
    end
    
    
    minImg = ((strongModel.meanImg - avgModel.minImg)*(j/numFrames)) + avgModel.minImg;
    varImg = ((strongModel.varImg - avgModel.varImg)*(j/numFrames)) +  avgModel.varImg;
     
    synFrame = normrnd(minImg,sqrt(varImg));
    synFrame = max(minImg, synFrame);
%     figure(i);
%     imshow(synFrame/255)

    outVariables.id = outVariables.id + 1;
    fprintf('Save Frame# %d / %d (%d)\n', ii, numFrames, outVariables.id);
    saveFrames(synFrame, outVariables)
end

end