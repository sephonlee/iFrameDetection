function finalFramesInfo = sampleFrames(iFramesInfo, streamLength, sectionSize, startShift)

importantPoint = iFramesInfo.importantPoint;
finalFramesInfo.sectionSize = sectionSize;
finalFramesInfo.startShift = startShift;
finalFramesInfo.streamLength = streamLength;
finalFramesInfo.finalFrames = [];
count = 0;

for i=1:1:size(importantPoint,1)
    
    if importantPoint(i, 1) - startShift < 1
        startPoint = 1;
    else
        startPoint = importantPoint(i, 1) - startShift;
    end
    
    if importantPoint(i, 1) + sectionSize - startShift - 1 > streamLength
        endPoint = streamLength;
    else
        endPoint = importantPoint(i, 1) + sectionSize - startShift - 1;
    end
    
    if isempty(finalFramesInfo.finalFrames)
        finalFramesInfo.finalFrames = [finalFramesInfo.finalFrames; [startPoint endPoint]];
        count = count + 1;
    else
        if finalFramesInfo.finalFrames(count, 2) >= startPoint;
            finalFramesInfo.finalFrames(count, 2) = endPoint;
        else
            finalFramesInfo.finalFrames = [finalFramesInfo.finalFrames; [startPoint endPoint]];
            count = count + 1;
        end
    end
end
end