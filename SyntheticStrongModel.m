clear all;
directory = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/frames';
outPath = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/synModel/';
% directory = '//psf/Home/Desktop/Research/PNNL/corpus/2012 03 04  - Copy2';

workingDir = '//psf/Home/Desktop/Research/PNNL/video/';
dirData = dir(directory);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
dirIndex(1,3) = 1;
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files



thresVar = 5;

height = 836;
width = 1232;
streamLength = length(fileList);

meanImg = zeros(height, width);
varImg = zeros(height,width);
minImg = zeros(height,width);
tempMat = [];
count = 0;
for w = 1:28:1232
    
    
    count = count + 1;
    tempMat = [];
    
    countImg = 0;
%     for ii = 279:291  strongest
        for ii = 174:179
        
        countImg = countImg +1;
        fprintf('Set #%d, frame #%d\n', count, ii);
       
        img = imread(fullfile(directory,fileList{ii}));
        imgOri = img;
        img = rgb2gray(img);
        img = double(img);
        % Filter legends
        img(255:310,1110:1230) = 0;
        img(200:230,1110:1230) = 0;
        img(250:275,1:125)= 0;
        img(725:830,260:335) = 0;
        img(790:825,1100:1230) = 0;
%         img_ = img(245:end, 243:end);
        
        tempMat(:,:,countImg) = img(:, w:w+27);
        
    end
    
    meanImg(:, w:w+27) = mean(tempMat, 3);
    minImg(:, w:w+27) = min(tempMat, [],3);
    maxImg(:, w:w+27) = max(tempMat, [],3);
    varImg(:, w:w+27) = var(tempMat, 0, 3);
    
end

maskActPixel = zeros(height, width);
maskActPixel(find(varImg > thresVar)) = 1;
% imshow(maskActPixel)


save(strcat(outPath,'moderateModel.mat'), 'meanImg', 'minImg', 'maxImg', 'varImg', 'maskActPixel');


meanImg
% save(strcat(outPath,'minImg.mat'), 'minImg')
% save(strcat(outPath,'varImg.mat'), 'varImg')

meanImgAll = load (strcat(outPath,'meanImg.mat'));
minImgAll = load (strcat(outPath,'minImg.mat'));
varImgAll = load (strcat(outPath,'varImg.mat'));

save(strcat(outPath,'avgModel.mat'), 'meanImg', 'minImg', 'varImg');

strongModel = load (strcat(outPath,'strongestModel.mat'));
avgModel = load (strcat(outPath,'avgModel.mat'));
moderateModel = load (strcat(outPath,'moderateModel.mat'));


makeStrongSignal(strongModel, avgModel, numFrames, needShaffle, randWindow)
makeStrongSignal(strongModel, avgModel, 10, 1, 2);
makeStrongSignal(moderateModel, avgModel, 10, 1, 2);

testMinImg = ((meanImg - minImgAll.minImg)*0.5) + minImgAll.minImg;
testVarImg = ((varImg - varImgAll.varImg)*0.5) +  varImgAll.varImg;

synFrame = normrnd(meanImg,sqrt(varImg));
synFrame = max(minImgAll.minImg, synFrame);
imshow(synFrame/255)

figure(2)
synFrame = normrnd(testMinImg,sqrt(testVarImg));
synFrame = max(testMinImg, synFrame);
imshow(synFrame/255)