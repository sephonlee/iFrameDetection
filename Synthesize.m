outVariables.outPath = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/synGrayFrames/';
outVariables.id = 0;
outVariables.fullDigits = 8;
modelInfo.sourcePath = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/synModel/';
modelInfo.numFrames = 1000;
modelInfo.models = {'average', 'moderate', 'strong'};
modelInfo.proModel = [75, 15, 10];
modelInfo.numModel = length(modelInfo.models);
modelInfo.maxNumInSection = 50;
modelInfo.modelIDArray = [];
modelInfo.modelNumFramesArray = [];


strongModel = load (strcat(modelInfo.sourcePath,'strongModel.mat'));
avgModel = load (strcat(modelInfo.sourcePath,'avgModel.mat'));
moderateModel = load (strcat(modelInfo.sourcePath,'moderateModel.mat'));

totalFrames = modelInfo.numFrames;

while outVariables.id < totalFrames
    
    numFramesInThisSection = randi([1 modelInfo.maxNumInSection]);
    lottery = randi([1 100]);
    
    if lottery - modelInfo.proModel(1) - modelInfo.proModel(2) > 0
        modelID = 3;
        makeAnomalyFrame(strongModel, avgModel, numFramesInThisSection, outVariables, 0);
        
    elseif lottery - modelInfo.proModel(1) > 0
        modelID = 2;
        makeAnomalyFrame(moderateModel, avgModel, numFramesInThisSection, outVariables, 0);
    else
        modelID = 1;
        makeSynFrame(avgModel, outVariables, numFramesInThisSection);
    end
    
    modelInfo.modelIDArray = [modelInfo.modelIDArray; modelID];
    modelInfo.modelNumFramesArray = [modelInfo.modelNumFramesArray; numFramesInThisSection];
    
    
    outVariables.id = outVariables.id + numFramesInThisSection;
    
end

modelInfo.numFrames = outVariables.id;

