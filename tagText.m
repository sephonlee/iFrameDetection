function InsertedImage = tagText(img, ii, clipID)

textInFrame = sprintf('Frame #%d, Clip #%d', ii, clipID);
TI = vision.TextInserter(textInFrame);
TI.Color = [1.0 1.0 0];
TI.FontSize = 20;
[height width] = size(img);
TI.Location = [floor(width/4), 1];
InsertedImage = step(TI, img);
% imshow(InsertedImage);

end