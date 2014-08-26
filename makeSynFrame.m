function synFrame = makeSynFrame(avgModel, outVariables, numFrames, needShaffle)

if nargin == 2
    numFrames = 1;
elseif nargin == 4 && needShaffle == 1
    fprintf('shaffle model...\n');
    avgModel.varImg = shaffleImg(avgModel.varImg, avgModel.maskActPixel, 5);
    fprintf('New model completed...\n');
    
end

startFrameId = outVariables.id + 1;
endFrameId = outVariables.id + numFrames;

for i = 1:1:numFrames
    synFrame = normrnd(avgModel.meanImg,sqrt(avgModel.varImg));
    synFrame = max(avgModel.minImg, synFrame);
%     figure(i);
%     imshow(synFrame/255);
    outVariables.id = outVariables.id + 1;
    fprintf('Save Frame# %d / %d (%d)\n', i, numFrames, outVariables.id);
    saveFrames(synFrame, outVariables)
end

end