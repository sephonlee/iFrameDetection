clear all;
directory = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/frames';

workingDir = '//psf/Home/Desktop/Research/PNNL/video/';
dirData = dir(directory);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
dirIndex(1,3) = 1;
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files


% outputVideo = VideoWriter(fullfile(workingDir,'Parbellan.avi'));
% outputVideo.FrameRate = 15;
% open(outputVideo);
% 
% outputVideo2 = VideoWriter(fullfile(workingDir,'mask.avi'));
% outputVideo2.FrameRate = 15;
% open(outputVideo2);

t=2; %bounday buffer for mask
% kernel = fspecial('gaussian', [5 5]); 
kernel = fspecial('average', [5 5]);
co_varOfPixelValue = 2;
co_varofSDist = 5;
poolSize = 5;
thresDif = 10;
thresSimilarity = 0.51;

features = [];
scores = [];
clipStart = 1;
clipID = 1;
lastHistogram = zeros(1,16);
maxnumRefDifPixels = 0;
% lastBackgroundValue = 92;

for ii = 1:length(fileList)
    
    
    
    fprintf('frame #%d\n', ii);
    img = imread(fullfile(directory,fileList{ii}));
    img = rgb2gray(img);
%     img = img(5:401,2:272);
    img = double(img);
    maxImg = max(max(img));
    minImg = min(min(img));

    imgUniform = conv2(img, kernel);
    imgUniform = img;
   
    if ii == 1
        refImgUniform = imgUniform;
        lastImgUniform = imgUniform;
    end
    
    [height width] = size(imgUniform);
    numPixels = height*width;
    
    % Differnce Mat
    difMat = abs(imgUniform - lastImgUniform);
    numDifPixels = length(find(difMat > thresDif));
    lastImgUniform = imgUniform;
    
    % Ref Difference Mat
    refDifMat = abs(imgUniform - refImgUniform);
    numRefDifPixels = length(find(refDifMat > thresDif));
    
    
    % Varience
    varPixelValues = var(imgUniform/mean(imgUniform));%%%%
    % Histogram
    histogram = hist(reshape(imgUniform, [1,numPixels]), 16)/(numPixels);
    histogramDist = sum(min([lastHistogram', histogram'],[], 2));
    lastHistogram = histogram;
    
    features = [features; [numDifPixels, numRefDifPixels, varPixelValues, histogramDist]];
    
    % Calculate overall score
%     score = co_varOfPixelValue*varPixelValues + co_varofHistDist*histogramDist;
    score = numDifPixels;
    
    if ii == 1
        tempScorePool = ones(poolSize, 1) *  score;
        refScore = score;
        % Accumalative Score Change
        lastScore = score;
        acumScoreChange = 0;   
    else
        tempScorePool = [tempScorePool(2:end); score];
    end
       
    acumScoreChange = acumScoreChange + abs(score - lastScore);
    lastScore = score;
       
    avgScore = mean(tempScorePool);
    
    magnitude = sum(sum(imgUniform))/numPixels;
    
    scores = [scores; [score acumScoreChange numRefDifPixels magnitude]];
    
    maxnumRefDifPixels = max(numRefDifPixels, maxnumRefDifPixels);
    
    if numRefDifPixels/numPixels > thresSimilarity
        clipStart = [clipStart; ii];
        refImgUniform = imgUniform;
        acumScoreChange = 0;
        clipID = clipID + 1;
    end
    
    img = img/255;
%     imgROI = img(1300:1800,1850:2250);
%     writeVideo(outputVideo,tagText(img, ii, clipID));
%     writeVideo(outputVideo2,mask);

end

% close(outputVideo)
% close(outputVideo2)

% figure(1)
% plot(1:length(fileList),features(:,1), '-ro',1:length(fileList), features(:,2),  '-rx', 1:length(fileList), features(:,3), '-b.', 1:length(fileList), features(:,4), '-b*', 1:length(fileList), features(:,5), '-gd');
% hleg1 = legend('Varience of PixelValues in ROI','Number of Pixels in ROI', 'Varience of Silhouette Distance', 'Average of Silhouette Distance', 'Distance of Histogram');

% [centroids class] = kmeans(features',3);
% hold on;
% plot(1:length(fileList),class(1,:)*100, 'c^')


figure(2)
plot(scores(:,1))
hleg1 = legend('Number of Difference Pixel', 'Accumalative Number of Difference Pixel', 'Number of Difference Pixel to Reference Frame');
for i=1:1:length(clipStart)
    hold on;
    plot([clipStart(i), clipStart(i)], [0, max(scores(:,2))], '--k')
end

[centroids class] = kmeans(scores(:,3)',3);
hold on;
plot(1:length(fileList),class(1,:)*10000, 'c^')
