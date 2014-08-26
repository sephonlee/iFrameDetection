
allResults = zeros(680,8);
for i=6:1:10

dir = sprintf('synShaffleFrames%d/', i);
dir = strcat('//psf/Home/Desktop/Research/PNNL/corpus/Nature/corpus/', dir);
results = dlmread(strcat(dir, 'allResult.txt'),'\t');

allResults = allResults+results;
end

allResults = allResults/ 5;

tailSizes = [30 50 80 100 200 500 1000 2000];
sectionSizes = [10 15 20 25 30 60];

for i = 1:1:length(sectionSizes)
   
    subList = allResults(find(allResults(:,2) == sectionSizes(i)), :);
    
    
    subList = allResults(find((allResults(:,2) == 25) & allResults(:,3) == 0), :);
end

subList2 = allResults(find((allResults(:,1) == 100) & (allResults(:,3) == 0)), :);

subList3 = allResults(find((allResults(:,1) == 100) & (allResults(:,2) == 30)), :);
