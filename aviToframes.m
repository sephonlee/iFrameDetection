outPath = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/framesGrey/';
file = '//psf/Home/Desktop/Research/PNNL/corpus/Nature/nmeth.2434-sv4.avi';

obj = mmreader(file);
vid = read(obj);
frames = obj.NumberOfFrames;
fullDigits = numel(num2str(frames));
for x = 1 : frames
    
    
    numAddDigits = fullDigits - numel(num2str(x));
    
    nameFile = num2str(x);
    if numAddDigits > 0
        for i=1:numAddDigits
            nameFile = strcat('0', nameFile);
        end
    end
    
   
    imwrite(rgb2gray(vid(:,:,:,x)),strcat(outPath, 'frame-', nameFile,'.jpg'));
end
