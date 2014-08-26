clear all;
% sumInfo.directory = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/synGrayFrames';
sumInfo.directory = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/synShaffleFrames5/';
modelInfo = load(strcat(sumInfo.directory,'modelInfo.mat'));
modelInfo = modelInfo.modelInfo;
frameDir = strcat(sumInfo.directory,'frames');
% sumInfo.directory = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/frames';
% sumInfo.directory = '//psf/Home/Desktop/Research/PNNL/corpus/2012 03 04  - Copy2';

workingDir = '//psf/Home/Desktop/Research/PNNL/video/';
dirData = dir(frameDir);      %# Get the data for the current sumInfo.directory
dirIndex = [dirData.isdir];  %# Find the index for directories
% dirIndex(1,3) = 1;
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files




% sumInfo.kernel = fspecial('gaussian', [5 5]);
sumInfo.kernel = fspecial('average', [5 5]);
% co_varOfPixelValue = 2;
% co_varofSDist = 5;
% poolSize = 5;
sumInfo.thresDif = 10;
% thresSimilarity = 0.51;

sumInfo.segments = {};
sumInfo.features = [];
sumInfo.segmentFeatures = [];
scores = [];
clipStart = 1;
clipID = 1;
sumInfo.clipDize = 50;
lastHistogram = zeros(1,16);
maxnumRefDifPixels = 0;
% lastBackgroundValue = 92;
sumInfo.streamLength = length(fileList);

sumInfo.streamLength = 5000;
count = 0;

for ii = 1:sumInfo.streamLength
    
    fprintf('frame #%d\n', ii);
    count = count + 1;
    img = imread(fullfile(frameDir,fileList{ii}));
    imgOri = img;
    img = rgb2gray(img);
    %     img = img(5:401,2:272);
    img = double(img);
    img = filterTimmer(img);
    maxImg = max(max(img));
    minImg = min(min(img));
    
    imgUniform = conv2(img, sumInfo.kernel);
    %     imgUniform = img;
    
    if ii == 1
        refImgUniform = imgUniform;
        lastImgUniform = imgUniform;
        avgImg = imgUniform;
        numDifPixelsArray = [];
        numRefDifPixelsArray = [];
    end
    
    [height width] = size(imgUniform);
    numPixels = height*width;
    
    % Differnce Mat
    difMat = abs(imgUniform - lastImgUniform);
    numDifPixels = length(find(difMat > sumInfo.thresDif));
    lastImgUniform = imgUniform; %update lastImg to current Img
    numDifPixelsArray = [numDifPixelsArray; numDifPixels];
    
    % Ref Difference Mat
    refDifMat = abs(imgUniform - refImgUniform);
    numRefDifPixels = length(find(refDifMat > sumInfo.thresDif));
    numRefDifPixelsArray = [numRefDifPixelsArray; numRefDifPixels];
    
    intensityGain = sum(sum(refDifMat)/height/width);
    
    % Varience
    varPixelValues = var(imgUniform/mean(imgUniform));%%%%
    % Histogram
    histogram = hist(reshape(imgUniform, [1,numPixels]), 16)/(numPixels);
    histogramDist = sum(min([lastHistogram', histogram'],[], 2));
    lastHistogram = histogram;
    
    
    if (mod(ii,sumInfo.clipDize) == 1) % start of a clip
        % Initialize a new Clip Info
        %         avgImg = imgUniform;
        greatestNumRefDifPixels = 0;
        characterImgUniform = imgUniform;
        characterImgID = ii;
        characterImg = imgOri;
    elseif (mod(count,sumInfo.clipDize) == 0) % end of a clip
        % Save Clip Info

        fprintf('Clip #%d\n', clipID);
        % Update refImageUniform for next clip 
        refImgUniform = avgImg;
%         figure(ii)
%         imshow(refImgUniform/255);
        
        % Save Clip Info 
        segmentInfo.startFrame = ii - sumInfo.clipDize + 1;
        segmentInfo.endFrame = ii;
        segmentInfo.avgImg = avgImg;
        segmentInfo.characterImgUniform = characterImgUniform;
        segmentInfo.characterImgID = characterImgID;
        segmentInfo.characterImg = characterImg;
        segmentInfo.difScore = greatestNumRefDifPixels;
        
       
        segmentInfo.numDifPixelsArray = numDifPixelsArray;
        segmentInfo.numRefDifPixelsArray = numRefDifPixelsArray;
        segmentInfo.varNumRefDifPixel = var(numRefDifPixelsArray);
        segmentInfo.varNumDifPixel = var(numDifPixelsArray);
        
        
        sumInfo.segments{clipID} = segmentInfo;
        clipID = clipID + 1;
        numDifPixelsArray = [];
        numRefDifPixelsArray = [];
        sumInfo.segmentFeatures = [sumInfo.segmentFeatures; [segmentInfo.difScore  segmentInfo.varNumRefDifPixel segmentInfo.varNumDifPixel]];
        
    else
        
        
%         avgImg = imgUniform.*(ii-1)/ii + avgImg./ii;
        
        % Update character of clip
        if numRefDifPixels > greatestNumRefDifPixels
            greatestNumRefDifPixels = numRefDifPixels;
            characterImgUniform = imgUniform;
            characterImg = imgOri;
            characterImgID = ii;
        end
         
    end
    
    avgImg = imgUniform.*(ii-1)/ii + avgImg./ii;
    sumInfo.features = [sumInfo.features; [numDifPixels, numRefDifPixels, varPixelValues, histogramDist intensityGain]];
    
    
end

% sumInfo.segmentInfo = segmentInfo;
save(strcat(sumInfo.directory, 'sumInfo.mat'), 'sumInfo');

sumInfo = load(strcat(sumInfo.directory, 'sumInfo.mat'));
sumInfo = sumInfo.sumInfo;

% groundTruth = [174:178, 279:285];
groundTruth = [279:285];
% groundTruth = [174:178, 195, 226, 279:291, 375, 432, 455, 471];
% iFramesInfo = detectChangePoint(sumInfo.features(:,2));



tailSizes = [30 50 80 100 200 500 1000 2000];
sectionSizes = [10 15 20 25 30 60];
% DIMS = [64 128 256];
% NC = [400,300];

allResults = [];
rfd = combvec(sectionSizes, tailSizes);

for i = 1:1:size(rfd, 2)
    
    sectionSize = rfd(1,i);
    tailSize = rfd(2,i);
    
    for startShift = 0:2:sectionSize
        
        fprintf('tailSize: %d, sectionSize: %d, startShift: %d\n', tailSize, sectionSize, startShift);
        iFramesInfo = detectChangePointBeta(sumInfo.features(:,2), tailSize);
        finalFramesInfo = sampleFrames(iFramesInfo, sumInfo.streamLength, sectionSize, startShift);
        evaResult =  evaluation(iFramesInfo, modelInfo, groundTruth, finalFramesInfo);
        
        
        allResults = [allResults; [tailSize, sectionSize, startShift, evaResult.recall(1), evaResult.precision(2), evaResult.compressionRate, evaResult.gtCompressionRate, iFramesInfo.maxSize]];        
    end
end

dlmwrite(strcat(sumInfo.directory, 'allResult.txt'), allResults, '\t');

%  iFramesInfo = detectChangePointBeta(sumInfo.features(:,2), 100, 1);
% finalFramesInfo = sampleFrames(iFramesInfo, sumInfo.streamLength, 20, 4);
%         evaResult =  evaluation(iFramesInfo, modelInfo, groundTruth, finalFramesInfo);
        

figure(2)
subplot(2,1,1)
hold on;
plot (evaResult.gtiFrames, sumInfo.features(evaResult.gtiFrames, 2), 'sc');
plot (evaResult.gtfFrames, sumInfo.features(evaResult.gtfFrames, 2), '*y');








% figure(1)
% plot(1:length(fileList),features(:,1), '-ro',1:length(fileList), features(:,2),  '-rx', 1:length(fileList), features(:,3), '-b.', 1:length(fileList), features(:,4), '-b*', 1:length(fileList), features(:,5), '-gd');
% hleg1 = legend('Varience of PixelValues in ROI','Number of Pixels in ROI', 'Varience of Silhouette Distance', 'Average of Silhouette Distance', 'Distance of Histogram');

% [centroids class] = kmeans(features',3);
% hold on;
% plot(1:length(fileList),class(1,:)*100, 'c^')
% 

charFrames = [];
clipStarts = [];
difScores = [];
for i = 1:1:length(sumInfo.segments)
    charFrames = [charFrames; sumInfo.segments{1,i}.characterImgID];
    difScores = [difScores; sumInfo.segments{1,i}.difScore];
    clipStarts = [clipStarts; sumInfo.segments{1,i}.startFrame];
end


Y = reshape(sumInfo.features(:,2),[sumInfo.clipDize, length(clipStarts)]);
Y = Y';

t = 1:sumInfo.clipDize;
X = repmat(t, length(clipStarts), 1);


Yps = [];
coefs = [];
figure(1)
plot(sumInfo.features(:,2))
legend('Number of Difference Pixel to Reference Frame');
for i=1:1:length(clipStarts)
    hold on;
    plot([clipStarts(i), clipStarts(i)], [0, max(difScores)], '--k')
    hold on;
    scatter(charFrames(i), sumInfo.features(charFrames(i),2), 'rx');
    
    [p,S] = polyfit(X(i,:),Y(i,:),4);
    %     p(5) = 0;
    Yp = polyval(p,X(i,:));
    
    Yps = [Yps; Yp];
    coefs = [coefs; p];
end

Yps = reshape(Yps', [sumInfo.streamLength,1]);
hold on ;
plot(Yps, 'r')

hold on;
plot(floor(sumInfo.clipDize/2):sumInfo.clipDize:sumInfo.streamLength, sumInfo.segmentFeatures(:,1), 'c^')





[centroids class] = kmeans(coefs(:,1:end)',5);
hold on;
plot(floor(sumInfo.clipDize/2)+100:sumInfo.clipDize:100+length(class)*sumInfo.clipDize,class(1,:)*20000, 'k^')



avgData = sumInfo.segments{1,1}.numRefDifPixelsArray;
m = ar(avgData, 25, 'yw');
allDif = [];
for i = 2:1:length(clipStarts)
    
    data = sumInfo.segments{1,i}.numRefDifPixelsArray;
    avgData = (i-1)*avgData./i + data./i;
    predData = predict(m, data, 1);
    m = ar(avgData, 25, 'yw');
    dif = sum(abs(data(1:49) - predData(2:50)));
    allDif = [allDif; dif];
end
    
m = ar(data, p, 'yw');    % yw for Yule-Walker method
pred = predict(m, data, 1);

coeffs = m.a;
nextValue = pred(end);

subplot(121), plot(data)
subplot(122), plot( (pred) )


figure(2)
plot(scores(:,1))


hleg1 = legend('Number of Difference Pixel', 'Accumalative Number of Difference Pixel', 'Number of Difference Pixel to Reference Frame');
for i=1:1:length(clipStart)
    hold on;
    plot([clipStart(i), clipStart(i)], [0, max(scores(:,2))], '--k')
end

[centroids class] = kmeans(scores(:,3)',3);
hold on;
plot(1:sumInfo.streamLength,class(1,:)*10000, 'c^')
