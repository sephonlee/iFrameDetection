function result = evaluation(iFramesInfo, modelInfo, groundTruth, finalFramesInfo)


% groundTruth = [174:178, 195, 226, 279:291, 375, 432, 455, 471];
% groundTruth = [174:178, 279:285];
groundTruth = sort(groundTruth);



aucList = modelInfo.modelFrameInfo(1:length(iFramesInfo.P)-1,:);
aucList(:,3:4) = 0;


for i = 1:1:size(aucList, 1)
    % is important point
    if sum(aucList(i, 2) == groundTruth)
        
        % mark as importand point
        aucList(i, 3) = 1;
    end
    
    index = 1;
    while (aucList(i,1) - finalFramesInfo.finalFrames(index, 2) > 0) && (index < size(finalFramesInfo.finalFrames, 1))
        index = index + 1;
    end
    
    if aucList(i,1) >= finalFramesInfo.finalFrames(index, 1) && aucList(i,1) <= finalFramesInfo.finalFrames(index, 2)
        % mark as true positive
        aucList(i, 4) = index;
    end

end



result.gTrueFrames = find(aucList(:, 3) == 1);
result.gFalseFrames = find(aucList(:, 3) == 0);
result.dTrueFrames = find(aucList(:, 4) > 0);
result.dFalseFrames = find(aucList(:, 4) == 0);
result.FNFrames = find((aucList(:, 3) == 1) & (aucList(:, 4) == 0));


correctTrue = length(find((aucList(:, 3) == 1) & (aucList(:,4) > 0)));
correctFalse = length(find((aucList(:, 3) == 0) & (aucList(:,4) == 0)));

result.recall(:,1) = correctTrue/length(result.gTrueFrames);
result.recall(:,2) = correctFalse/length(result.gFalseFrames);


result.precision(:,1) = correctTrue/length(result.dTrueFrames);

if correctFalse == 0
    result.precision(:,2) = 1;
else
    result.precision(:,2) = correctFalse/length(result.dFalseFrames);
end

totalFinalFrame = 0;
for i=1:1:size(finalFramesInfo.finalFrames,1)
    totalFinalFrame = totalFinalFrame + (finalFramesInfo.finalFrames(i, 2) - finalFramesInfo.finalFrames(i, 1) + 1);
end

result.compressionRate = totalFinalFrame/size(modelInfo.modelFrameInfo,1);


iFramesTruth.importantPoint = result.gTrueFrames;

finalTruthFramesInfo = sampleFrames(iFramesTruth, finalFramesInfo.streamLength, finalFramesInfo.sectionSize, finalFramesInfo.startShift);

totalFinalFrame = 0;
for i=1:1:size(finalTruthFramesInfo.finalFrames,1)
    totalFinalFrame = totalFinalFrame + (finalTruthFramesInfo.finalFrames(i, 2) - finalTruthFramesInfo.finalFrames(i, 1) + 1);
end

result.gtCompressionRate = totalFinalFrame/size(modelInfo.modelFrameInfo,1);


end