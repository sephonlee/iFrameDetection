
directory = '//psf/Home/Desktop/Research/PNNL/corpus/Ideal/';
workingDir = '//psf/Home/Desktop/Research/PNNL/video/';
dirData = dir(directory);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
dirIndex(1,3) = 1;
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files


% outputVideo = VideoWriter(fullfile(workingDir,'partial.avi'));
% outputVideo.FrameRate = 15;
% open(outputVideo);
% 
% outputVideo2 = VideoWriter(fullfile(workingDir,'mask.avi'));
% outputVideo2.FrameRate = 15;
% open(outputVideo2);

t=2; %bounday buffer for mask
co_varOfPixelValue = 2;
co_varofSDist = 5;
poolSize = 5;
thresSimilarity = 0.15;

features = [];
scores = [];
clipStart = 1;
clipID = 1;
lastHistogram = zeros(1,16);

% lastBackgroundValue = 92;

for ii = 1:length(fileList)
    fprintf('frame #%d\n', ii);
    img = imread(fullfile(directory,fileList{ii}));
    img = double(img);
    maxImg = max(max(img));
    minImg = min(min(img));
    
%     img = img/maxImg;
    
    % Zoom into ROI
    varH = var(img, [], 1);
    targetIndice = find(varH > 10e-5);
    minCol = min(targetIndice);
    maxCol = max(targetIndice);
    
    varV = var(img, [], 2);
    targetIndice = find(varV > 10e-5);
    minRow = min(targetIndice);
    maxRow = max(targetIndice);
    
    % normalized illumination
    backgroundValue = img(floor(minRow/2), floor(minCol/2));
    if ii == 1
        lastBackgroundValue = backgroundValue;
    end 
    img = img/backgroundValue*lastBackgroundValue;
    img(find(img > 255)) = 255;
    backgroundValue = img(floor(minRow/2), floor(minCol/2));
    % find ROI
    mask = zeros(size(img));
    mask(find(img > backgroundValue + t | img < backgroundValue - t)) = 1;
    % filter noise
    se = strel('disk',3);
    mask = imclose(mask, se);
    
    [row col] = find(mask == 1);
    pixelValuesOfROI = img(find(mask == 1));
    varPixelValuesOfROI = var(pixelValuesOfROI/mean(pixelValuesOfROI));%%%%
    numPixelValuesOfROIper100 = length(row)/100;
    
    % find contour of ROI
    contour = bwtraceboundary(mask, [floor((minRow+maxRow)/2), floor((minCol+maxCol)/2)], 'W', 8, 50,...
        'counterclockwise');
    B = bwboundaries(mask);
    point = B{1};
    mass = mean(point, 1);
    silhouetteDist = sqrt(sum(bsxfun(@minus, mass, point).^2, 2));
    silhouetteDistAvg = mean(silhouetteDist);%%%%
    silhouetteDistVar = var(silhouetteDist/silhouetteDistAvg);%%%%

    
    % Histogram
    histogram = hist(pixelValuesOfROI, 16)/length(row);
    histogramDist = sum(min([lastHistogram', histogram'],[], 2));
    lastHistogram = histogram;
    
    features = [features; [varPixelValuesOfROI, numPixelValuesOfROIper100, silhouetteDistVar, silhouetteDistAvg, histogramDist*100]];
    
    % Calculate overall score
    co_varOfPixelValue = 2;
    co_varofSDist = 5;
    score = co_varOfPixelValue*varPixelValuesOfROI + co_varofSDist*silhouetteDistVar;
%     scores = [scores; [varPixelValuesOfROI silhouetteDistVar score]];
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
    scores = [scores; [varPixelValuesOfROI silhouetteDistVar score acumScoreChange]];
    
    if abs(avgScore - refScore)/refScore > thresSimilarity
        clipStart = [clipStart; ii];
        refScore = avgScore;
        acumScoreChange = 0;
        clipID = clipID + 1;
    end
    
    img = img/255;
    imgROI = img(1300:1800,1850:2250);
%     writeVideo(outputVideo,tagText(imgROI, ii, clipID));
%     writeVideo(outputVideo2,mask);

end

% close(outputVideo)
% close(outputVideo2)

figure(1)
plot(1:length(fileList),features(:,1), '-ro',1:length(fileList), features(:,2),  '-rx', 1:length(fileList), features(:,3), '-b.', 1:length(fileList), features(:,4), '-b*', 1:length(fileList), features(:,5), '-gd');
hleg1 = legend('Varience of PixelValues in ROI','Number of Pixels in ROI', 'Varience of Silhouette Distance', 'Average of Silhouette Distance', 'Distance of Histogram');

[centroids class] = kmeans(features',3);
hold on;
plot(1:length(fileList),class(1,:)*100, 'c^')


figure(2)
plot(scores)
hleg1 = legend('Varience of PixelValues in ROI', 'Varience of Silhouette Distance', 'Overall Score', 'Accumlative Score Change');
hold on
for i=1:1:length(clipStart)

    plot([clipStart(i), clipStart(i)], [0, max(scores(:,4))], '--k')
end



textInFrame = sprintf('Frame #%d, Clip #%d', ii, clipID);
TI = vision.TextInserter(textInFrame);
TI.Color = [1.0 1.0 0];
TI.FontSize = 20;
[height width] = size(imgROI);
TI.Location = [floor(width/4), 1];
InsertedImage = step(TI, imgROI);
imshow(InsertedImage);



% Zoom into ROI
varH = var(img, [], 1);
targetIndice = find(varH > 10e-5);
minCol = min(targetIndice);
maxCol = max(targetIndice);

varV = var(img, [], 2);
targetIndice = find(varV > 10e-5);
minRow = min(targetIndice);
maxRow = max(targetIndice);

% find ROI
backgroundValue = img(floor(minRow/2), floor(minCol/2));
mask = zeros(size(img));
mask(find(img ~= backgroundValue)) = 1;
% filter noise
se = strel('disk',3);
mask = imclose(mask, se);

[row col] = find(mask == 1);
pixelValuesOfROI = img(find(mask == 1));
varPixelValuesOfROI = var(pixelValuesOfROI);%%%%

% find contour of ROI
contour = bwtraceboundary(mask, [floor((minRow+maxRow)/2), floor((minCol+maxCol)/2)], 'W', 8, 50,...
    'counterclockwise');
B = bwboundaries(mask);
point = B{1};
mass = mean(point, 1);
silhouetteDist = sqrt(sum(bsxfun(@minus, mass, point).^2, 2));
silhouetteDistVar = var(silhouetteDist);%%%%
silhouetteDistAvg = mean(silhouetteDist);%%%%

features = [features; [varPixelValuesOfROI, silhouetteDistVar, silhouetteDistAvg]];

% imshow(mask)
% hold on;
% point = B{1};
% cVec1 = 'bgrcmykbgrcmykbgrcmykbgrcmykmykbgrcmykbgrcmykrcmykbgrcmykmykbgrcmykgrcmykbgrcmykbgrcmykmykbgrcmykbgrcmykrcmykbgrcmykmykbgrcmy';
% cVec2 = 'od+*^v<>phxsod+*^*^v<>phxsod+*^v<>phxsod+*^*^v<>phxsod+*^v<>phxsod+*^*^v<>phxd+*^*^v<>phxsod+*^v<>phxsod+*^*^v<sod+*^v<>phxsod';
% plot(point(:,2),point(:,1),'.','MarkerEdgeColor',cVec1(1));







imgR = img(minRow:maxRow,minCol:maxCol);