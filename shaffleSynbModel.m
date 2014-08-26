clear all;
outVariables.outPath = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/synShaffleFrames10/';
outVariables.id = 0;
outVariables.fullDigits = 8;

modelInfo.sourcePath = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/frames/';
modelInfo.numFrames = 8000;
modelInfo.minNumInSection = 80;
modelInfo.maxNumInSection = 120;
modelInfo.modelSectionInfoColumn = {'Start Frame ID', 'Number of Frames', 'Direction'};
modelInfo.modelSectionInfo = [];
modelInfo.modelFrameInfo = [];

mkdir(outVariables.outPath);
mkdir(strcat(outVariables.outPath, 'frames'));

dirData = dir(modelInfo.sourcePath);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
dirIndex(1,3) = 1;
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files

modelInfo.sourceStreamLength = length(fileList);
% modelInfo.sourceStreamLength = 500;
modelInfo.startSourceFrame = 1;

totalFrames = modelInfo.numFrames;

while outVariables.id < totalFrames
    
%     while numFrames >150 || numFrames <400
    sampleFrameID = randi([modelInfo.startSourceFrame modelInfo.sourceStreamLength]);
%     end
%     sampleFrameID = randi([1 200]);
    
    numFrames = randi([modelInfo.minNumInSection modelInfo.maxNumInSection]);
    direction = randi([1 2]);


    startFrame = sampleFrameID;
    if direction == 1
        endFrame = startFrame + numFrames - 1;
    else
        endFrame = startFrame - numFrames + 1;
    end
    
    if endFrame > modelInfo.sourceStreamLength
        endFrame = modelInfo.sourceStreamLength;
        %%%
        while endFrame - startFrame + 1 < modelInfo.minNumInSection
            fprintf('redo, %d\n', endFrame - startFrame + 1);
            sampleFrameID = randi([modelInfo.startSourceFrame modelInfo.sourceStreamLength]);
            startFrame = sampleFrameID;
            endFrame = startFrame + numFrames - 1;
            
            if endFrame > modelInfo.sourceStreamLength
                endFrame = modelInfo.sourceStreamLength;
            end
        end    
        %%%   
        numFrames = abs(startFrame - endFrame) + 1;
    elseif endFrame < modelInfo.startSourceFrame
        endFrame = modelInfo.startSourceFrame;
        %%%
        while startFrame - endFrame + 1 < modelInfo.minNumInSection
            fprintf('redo, %d\n', startFrame - endFrame + 1);
            sampleFrameID = randi([modelInfo.startSourceFrame modelInfo.sourceStreamLength]);
            startFrame = sampleFrameID;
            endFrame = startFrame - numFrames + 1;
            
            if endFrame < modelInfo.startSourceFrame;
                endFrame = modelInfo.sourceStreamLength;
            end
        end   
        %%%
 
        numFrames = abs(startFrame - endFrame) + 1;
    end
    
    modelInfo.modelSectionInfo = [modelInfo.modelSectionInfo; [outVariables.id + 1, sampleFrameID, numFrames, direction]];
    
    
    
    for i = 0:1:numFrames - 1
        
        if startFrame < endFrame
            ii = startFrame + i;
        else
            ii = startFrame - i;
        end
        
        img = imread(fullfile(modelInfo.sourcePath, fileList{ii}));
        
        outVariables.id = outVariables.id + 1;
        numAddDigits = outVariables.fullDigits - numel(num2str(outVariables.id));
        
        % Save frames
        nameFile = num2str(outVariables.id);
        if numAddDigits > 0
            for x = 1:numAddDigits
                nameFile = strcat('0', nameFile);
            end
        end
        
        % Record Frame ID
        modelInfo.modelFrameInfo = [modelInfo.modelFrameInfo; [outVariables.id ii]];

        fprintf('Save Frame# %d / %d, id: #%d (%d)\n', i + 1, numFrames, outVariables.id, ii);
        imwrite(img, strcat(outVariables.outPath, '/frames/frame-', nameFile,'.jpg'));
    end
    
end

save(strcat(outVariables.outPath,'modelInfo.mat'), 'modelInfo')