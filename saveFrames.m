function saveFrames(img, outVariables)

numAddDigits = outVariables.fullDigits - numel(num2str(outVariables.id));

nameFile = num2str(outVariables.id);
if numAddDigits > 0
    for i=1:numAddDigits
        nameFile = strcat('0', nameFile);
    end
end


imwrite(img/255, strcat(outVariables.outPath, 'frame-', nameFile,'.jpg'));

end