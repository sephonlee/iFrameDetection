clear all;
directory = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/frames';
% directory = '//psf/Home/Desktop/Research/PNNL/corpus/2012 03 04  - Copy2';

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

segments = {};
features = [];
scores = [];
clipStart = 1;
clipID = 1;
clipSize = 50;
lastHistogram = zeros(1,16);
maxnumRefDifPixels = 0;
% lastBackgroundValue = 92;
count = 0;
streamLength = length(fileList);
% streamLength=5000;

for ii = 1:streamLength
    
    fprintf('frame #%d\n', ii);
    count = count + 1;
    img = imread(fullfile(directory,fileList{ii}));
    imgOri = img;
    img = rgb2gray(img);
    img = img(5:401,2:272);
    img = double(img);
    maxImg = max(max(img));
    minImg = min(min(img));
    
    imgUniform = conv2(img, kernel);
    %     imgUniform = img;
    
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
    
    if (count == 1)
        avgImg = imgUniform;
        greatestNumRefDifPixels = 0;
        characterImgUniform = imgUniform;
        characterImgID = ii; 
        characterImg = imgOri;
    else
        avgImg = imgUniform.*(count-1)/count + avgImg./count;
        
        if numRefDifPixels > greatestNumRefDifPixels
            greatestNumRefDifPixels = numRefDifPixels;
            characterImgUniform = imgUniform;
            characterImg = imgOri;
            characterImgID = ii;
        end
        
        if(count == clipSize)
            fprintf('Clip #%d\n', clipID);
            % Update
            refImgUniform = avgImg;
            count = 0;
            
            % Save info
            difMat = abs(characterImgUniform - avgImg);
            numDifPixels = length(find(difMat > thresDif));
            
            segmentInfo.startFrame = ii - clipSize + 1;
            segmentInfo.endFrame = ii;
            segmentInfo.avgImg = avgImg;
            segmentInfo.characterImgUniform = characterImgUniform;
            segmentInfo.characterImgID = characterImgID;
            segmentInfo.characterImg = characterImg;
            segmentInfo.difScore = numDifPixels;
            
            segments{clipID} = segmentInfo;
            clipID = clipID + 1;
            
        end
    end
    
    
    features = [features; [numDifPixels, numRefDifPixels, varPixelValues, histogramDist]];
    
    
    % Calculate overall score
    %     score = co_varOfPixelValue*varPixelValues + co_varofHistDist*histogramDist;
    %     score = numDifPixels;
    %
    %     if ii == 1
    %         tempScorePool = ones(poolSize, 1) *  score;
    %         refScore = score;
    %         % Accumalative Score Change
    %         lastScore = score;
    %         acumScoreChange = 0;
    %     else
    %         tempScorePool = [tempScorePool(2:end); score];
    %     end
    %
    %     acumScoreChange = acumScoreChange + abs(score - lastScore);
    %     lastScore = score;
    %
    %     avgScore = mean(tempScorePool);
    %
    % %     magnitude = sum(sum(imgUniform))/numPixels;
    %
    %     scores = [scores; [score acumScoreChange numRefDifPixels]];
    %
    %     maxnumRefDifPixels = max(numRefDifPixels, maxnumRefDifPixels);
    %
    %     if numRefDifPixels/numPixels > thresSimilarity
    %         clipStart = [clipStart; ii];
    %         refImgUniform = imgUniform;
    %         acumScoreChange = 0;
    %         clipID = clipID + 1;
    %     end
    %
    %     img = img/255;
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


charFrames = [];
clipStarts = [];
difScores = [];
for i = 1:1:length(segments)
    charFrames = [charFrames; segments{1,i}.characterImgID];
    difScores = [difScores; segments{1,i}.difScore];
    clipStarts = [clipStarts; segments{1,i}.startFrame];
end


Y = reshape(features(:,2),[clipSize, length(clipStarts)]);
Y = Y';

t = 1:clipSize;
X = repmat(t, length(clipStarts), 1);


Yps = [];
coefs = [];
figure(1)
plot(features(:,2))
legend('Number of Difference Pixel to Reference Frame');
for i=1:1:length(clipStarts)
    hold on;
    plot([clipStarts(i), clipStarts(i)], [0, difScores(i)], '--k')
    hold on;
    scatter(charFrames(i), features(charFrames(i),2), 'rx');
    
    [p,S] = polyfit(X(i,:),Y(i,:),4);
%     p(5) = 0;
    Yp = polyval(p,X(i,:));
    
    Yps = [Yps; Yp];
    coefs = [coefs; p];
end

Yps = reshape(Yps', [streamLength,1]);
hold on ;
plot(Yps, 'r')

[centroids class] = kmeans(coefs(11:20,1:end-1)',5);
hold on;
plot(floor(clipSize/2)+100:clipSize:100+length(class)*clipSize,class(1,:)*20000, 'k^')


% 
% m = ar(data, p, 'yw');    % yw for Yule-Walker method
% pred = predict(m, data, 1);
% 
% coeffs = m.a;
% nextValue = pred(end);
% 
% subplot(121), plot(data)
% subplot(122), plot( cell2mat(pred) )


figure(2)
plot(scores(:,1))


hleg1 = legend('Number of Difference Pixel', 'Accumalative Number of Difference Pixel', 'Number of Difference Pixel to Reference Frame');
for i=1:1:length(clipStart)
    hold on;
    plot([clipStart(i), clipStart(i)], [0, max(scores(:,2))], '--k')
end

[centroids class] = kmeans(scores(:,3)',3);
hold on;
plot(1:streamLength,class(1,:)*10000, 'c^')
