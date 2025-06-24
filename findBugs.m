% for Azhar to identify the bugs in her image via DAPI and then use this to
% count how many red spots (a protein?) there are in the image

% outputs a table with columns
% cell number | cell centre x | cell centre y | number red spots |...
% Pearson correlation coefficient | spot1 x | spot 1 y | spot2 x | spot2 y
% ...

function findBugs()

% find the file we want to look at
[fileName,pathName] = uigetfile('*.tif');
image = [pathName,fileName];

% read in the image and select the DAPI image which is blue channel second
% image
T = tiffread(image);
filename = image(1:end-4);

DAPI = T(2).data;
DAPI = DAPI{3};

% find the bugs

% show the original image
disp(['max: ',num2str(max(DAPI(:))),' min: ',num2str(min(DAPI(:)))]);
figure;imshow(DAPI,[min(DAPI(:)) max(DAPI(:))]);title('original')

original_im = DAPI; % save the original image for later

min_px_in_bug = 100; %The minimum number of pixels per detected cell, which
%should remove the majority of background noise

max_px_in_bug = 2000; %The maximum area a detected cell can take, it
%should remove any large "blobs" of dirt on the sample

%Filtering the image using imageFilter and outputting the image to
%"filtered_im", this allows us to find the bugs
[filtered_im] = imageFilterRT(DAPI);

% find connected components
clusters = bwconncomp(filtered_im);

% remove clusters that are too big or too small to be bugs
numPixels = cellfun(@numel,clusters.PixelIdxList);

for i=1:length(clusters.PixelIdxList)
    if numPixels(i) > max_px_in_bug || numPixels(i) < min_px_in_bug
        filtered_im(clusters.PixelIdxList{i})=0;
    end
end

% find clusters again now we may have removed some

clusters = bwconncomp(filtered_im);

%figure;imshow(filtered_im,[min(filtered_im(:)),max(filtered_im(:))])

%%%%%%%

% find details about each of these clusters
s = regionprops(clusters,'BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Area','PixelList','Extrema','Centroid');
numBugs = numel(s);


% save these details into a table. This currently includes non-bugs and
% pairs of bugs
[bug_info_table] = zeros(numBugs,11);

for i=1:numBugs
    bug_info_table(i,1) = i; % first entry is bug id number
    bug_info_table(i,2) = s(i).Extrema(2,2); % top 
    bug_info_table(i,3) = s(i).Extrema(6,2); % bottom
    bug_info_table(i,4) = s(i).Extrema(2,1); % right
    bug_info_table(i,5) = s(i).Extrema(6,1); % left
    bug_info_table(i,6) = s(i).MajorAxisLength; % length long axis bug
    bug_info_table(i,7) = s(i).MinorAxisLength; % length short axis bug
    bug_info_table(i,6) = s(i).MajorAxisLength; % length long axis bug
    bug_info_table(i,8) = s(i).Orientation; % angle between x axis and long axis of bug
    bug_info_table(i,9) = i; % old bug ID number (currently same)
    
    % calculate the total intensity in the bug
    totalIntensity = 0;
    for j=1:numel(clusters.PixelIdxList{i})
        totalIntensity = totalIntensity + sum(original_im(clusters.PixelIdxList{i}(j)));
    end
    
    bug_info_table(i,10) = totalIntensity;

end

% display the original image
figure;imshow(original_im,[min(DAPI(:)) max(DAPI(:))]);

% and display a number next to each bug and the bounding box of the
% cluster that corresponds to that bug
for k = 1:numBugs
    hold on;

    rectangle('Position', s(k).BoundingBox,'EdgeColor','r', 'LineWidth', 1)
    text(s(k).BoundingBox(1)-5,s(k).BoundingBox(2)-5,num2str(k),'Color','g');
    %disp(['Num bug: ',num2str(i),' Long axis: ',num2str(s(i).MajorAxisLength),' Short axis: ',num2str(s(i).MinorAxisLength)])
end

% set all bugs to not included
bug_info_table(:,11) = 0; %None are included initially
numBugs = length(bug_info_table(:,1));

% pop up a checkbox where the bugs we want to analyse can be selected
chosenBugs = checkboxes(numBugs);

% close the figure so it's not cluttering our screen
close(gcf)

% now we set the bugs we want to analyse to 1
bug_info_table(:,11) = chosenBugs;

% get rid of all bugs that weren't chosen
bug_info_table(find(bug_info_table(:,11)==0),:) = [];

% find out how many bugs we need to analyse
numAnalysedBugs = size(bug_info_table,1);

% For ease of code save the old cluster number of each bug so we can refer
% back to its previous details
oldIndices = zeros(1,numAnalysedBugs);
for i=1:numAnalysedBugs
    oldIndices(i) = bug_info_table(i,1);
end

% and then give each bug a new number starting from 1
bug_info_table(:,1) = 1:size(bug_info_table,1);

% I now use the 11th column to record the old index of this bug
bug_info_table(:,11) = oldIndices;

% show the original image
figure;imshow(original_im,[min(DAPI(:)) max(DAPI(:))]);

% and show on top the bounding box and new number for each analysed bug
for k = 1:numAnalysedBugs
    hold on;
    oldK = bug_info_table(k,11);
    rectangle('Position', s(oldK).BoundingBox,'EdgeColor','r', 'LineWidth', 1)
    text(s(oldK).BoundingBox(1)-5,s(oldK).BoundingBox(2)-5,num2str(k),'Color','g');
    %disp(['Num bug: ',num2str(i),' Long axis: ',num2str(s(i).MajorAxisLength),' Short axis: ',num2str(s(i).MinorAxisLength)])
end

% save this figure
saveas(gcf,strcat(filename,'-analysedBugs.png'))

%%%%%

% now we need to find the red spots on each bug
% read in the red channel of the red image
red = T(3).data;
red = red{1};

% loop through each cell

% large data storage
spotData = zeros(numAnalysedBugs,30); % holds spot maxima x and y locations
spotIntensity = zeros(numAnalysedBugs,15); % holds the intensity of each spot
numSpots = zeros(numAnalysedBugs,1); % number of spots per cell
cellNumbers = (1:numAnalysedBugs)';
cellCentresX = zeros(numAnalysedBugs,1);
cellCentresY = zeros(numAnalysedBugs,1);
cellRadii = zeros(numAnalysedBugs,1);
PCC = zeros(numAnalysedBugs,1);
distances = zeros(numAnalysedBugs,15);
normDistances = zeros(numAnalysedBugs,15);

for i=1:numAnalysedBugs
    % make a mask of that cell
    % make a copy of filtered_im
    
    filtered_cell = zeros(size(filtered_im));
    oldK = bug_info_table(i,11);
    
    filtered_cell(clusters.PixelIdxList{oldK})=1;
    
    %figure;imshow(filtered_cell)
    
    % apply that to the red image
    
    redCell = red;
    redCell(filtered_cell==0)=0;
    
%     figure;imshow(redCell,[min(red(:)),max(red(:))])

    % compute the PCC
    redPixelVals = redCell(filtered_cell==1);
    bluePixelVals = DAPI(filtered_cell==1);
    
    PCC(i) = corr(double(redPixelVals(:)),double(bluePixelVals(:))) ;
    
    % median filter
    
    redCell = medfilt2(redCell,[2 2]);
    
    %figure;imshow(redCell,[min(red(:)),max(red(:))])
    
    % find the maxima in red image
    
    BW = imregionalmax(redCell);
    
    % find connected components, these are the peaks
    % connected components because some peaks might be more than one pixel
    maxima = bwconncomp(BW);
    
    maximaDetails = regionprops(maxima,'Centroid');
    centroids = cat(1, maximaDetails.Centroid);
    numMaxima = numel(maximaDetails);
    
    % check the size and replace large ones by the centre of mass
    
    onePxBW = zeros(size(BW));

    for k=1:numMaxima
        % replace each maxima by the centre of mass of each maxima
        % if the maxima is only one pixel it stays the same
        % if it's multiple it becomes its centre
        
        % the order of these coordinates is IMPORTANT
        onePxBW(floor(centroids(k,2)),floor(centroids(k,1))) = 1;
    end
    
%     figure;imshow(BW)
%     hold on;plot(centroids(:,1),centroids(:,2),'ro')
%     figure;imshow(onePxBW)
    
    linIndexes = find(onePxBW==1);
    
    % convert indexes to subscripts
    
    %[x,y] = ind2sub(size(red),linIndexes);
    
    %figure;imshow(redCell,[min(red(:)),max(red(:))])
    %hold on;
    %plot(y,x,'ro','MarkerSize',0.5)
    
    % show them to me
    % remove peaks that are below the average
    
    n = sum(redCell(:)~=0);
    m = sum(redCell(:)) ./ n;

    thresh = m;
    
    redCell(linIndexes)
    
    chosenPeaksLin = linIndexes(redCell(linIndexes)>thresh);
    
    [x,y] = ind2sub(size(red),chosenPeaksLin);
    
    %hold on;
    %plot(y,x,'go','MarkerSize',0.5)
    
    oldI = bug_info_table(i,11);
    cellCentre = s(oldI).Centroid; % 1st index is row index 2nd is column
    
    % save this information

    cellCentresX(i) = cellCentre(2);
    cellCentresY(i) = cellCentre(1);
    
    % old code where radius was average of major and minor axes
%     cellDiameter = 0.5*(s(oldI).MajorAxisLength+s(oldI).MinorAxisLength);
%     cellRadius = 0.5*cellDiameter;
    
    % the cell radius will be the maximum distance from centre
    % first create linear indices for each pixel in bug
    
    linIndices = clusters.PixelIdxList{oldK};
    % convert to subscript
    [bugX,bugY] = ind2sub(size(red),linIndices);
    
    maxDist = 0;
    maxX = 0;
    maxY = 0;
    
    for k=1:length(linIndices)
        centreDist = (bugX(k) - cellCentre(2))^2+(bugY(k) - cellCentre(1))^2;
        if centreDist > maxDist
            maxDist = centreDist;
            maxX = bugX(k);
            maxY = bugY(k);
        end
    end
    
    %plot(cellCentre(1),cellCentre(2),'bo')
    %plot(maxY,maxX,'yo')
    
    cellRadius = maxDist^0.5;
    cellRadii(i) = cellRadius;
    
    
    %plot(y(1),x(1),'bo','MarkerSize',0.5)
    
    %plot(floor(cellCentre(1)+0.5*(y(1)-cellCentre(1))),floor(cellCentre(2)+0.5*(x(1)-cellCentre(2))),'yo')
    
    % output
    
    % save the x and y coordinates of each maxima
    
    numPeaks = numel(x);
    
    for n=1:numPeaks
        spotData(i,(2*n-1)) = x(n);
        spotData(i,(2*n)) = y(n);
    end

    % count the number of spots
    
    numSpots(i) = numPeaks;
    
    % for each spot calculate the distance from the centre of the cell
    
    % cell centre was calculated above
    
    % y distance is y(n) - cell centre(1), x distance is x(n) - cell
    % centre(2)
    
    for n=1:numPeaks
        xDist = x(n)-cellCentre(2);
        yDist = y(n)-cellCentre(1);
        
        dist = sqrt(xDist^2+yDist^2);
        
        normDist = dist/cellRadius;
        
        distances(i,n) = dist;
        normDistances(i,n) = normDist;
        spotIntensity(i,n) = redCell(chosenPeaksLin(n));
    end

end

% plot all the peaks for reference

figure;imshow(red,[min(red(:)),max(red(:))])
hold on;
% add the peaks

for i=1:(size(spotData,1))
    for e=1:(size(spotData,2)/2)
        plot(spotData(i,2*e),spotData(i,2*e-1),'ro','MarkerSize',0.5)
    end
end

% save this figure
saveas(gcf,strcat(filename,'-foundSpots.png'))
close(gcf)

% warn if too many spots found on one bug
if size(spotData,2) > 30
    disp('WARNING: for one cell more than 15 spots were found, data not recorded properly')
    disp(['Max number of spots found was ',num2str(size(spotData,2)/2)])
end

% disp('spotdata')
% spotData  % holds spot maxima x and y locations
% disp('numspots')
% numSpots % number of spots per cell
% disp('cellnumbers')
% cellNumbers
% disp('cellcentresx')
% cellCentresX
% disp('cellcentresy')
% cellCentresY
% disp('pcc')
% PCC
% disp('distances')
% distances
% disp('normdistances')
% normDistances 
    
outputTable = table(cellNumbers,cellCentresX,cellCentresY,cellRadii,PCC,numSpots,...
    spotData(:,1),spotData(:,2),spotIntensity(:,1),distances(:,1),normDistances(:,1),...
    spotData(:,3),spotData(:,4),spotIntensity(:,2),distances(:,2),normDistances(:,2),...
    spotData(:,5),spotData(:,6),spotIntensity(:,3),distances(:,3),normDistances(:,3),...
    spotData(:,7),spotData(:,8),spotIntensity(:,4),distances(:,4),normDistances(:,4),...
    spotData(:,9),spotData(:,10),spotIntensity(:,5),distances(:,5),normDistances(:,5),...
    spotData(:,11),spotData(:,12),spotIntensity(:,6),distances(:,6),normDistances(:,6),...
    spotData(:,13),spotData(:,14),spotIntensity(:,7),distances(:,7),normDistances(:,7),...
    spotData(:,15),spotData(:,16),spotIntensity(:,8),distances(:,8),normDistances(:,8),...
    spotData(:,17),spotData(:,18),spotIntensity(:,9),distances(:,9),normDistances(:,9),...
    spotData(:,19),spotData(:,20),spotIntensity(:,10),distances(:,10),normDistances(:,10),...
    spotData(:,21),spotData(:,22),spotIntensity(:,11),distances(:,11),normDistances(:,11),...
    spotData(:,23),spotData(:,24),spotIntensity(:,12),distances(:,12),normDistances(:,12),...
    spotData(:,25),spotData(:,26),spotIntensity(:,13),distances(:,13),normDistances(:,13),...
    spotData(:,27),spotData(:,28),spotIntensity(:,14),distances(:,14),normDistances(:,14),...
    spotData(:,29),spotData(:,30),spotIntensity(:,15),distances(:,15),normDistances(:,15));

% change headings

outputTable.Properties.VariableNames = {'CellNumber','CellCentreX','CellCentreY','CellRadius'...
    'PCC','NumRedSpots',...
    'Spot1X','Spot1Y','Spot1Intensity','Spot1toCentreDist','Spot1toCentreDistNormalised',...
    'Spot2X','Spot2Y','Spot2Intensity','Spot2toCentreDist','Spot2toCentreDistNormalised',...
    'Spot3X','Spot3Y','Spot3Intensity','Spot3toCentreDist','Spot3toCentreDistNormalised',...
    'Spot4X','Spot4Y','Spot4Intensity','Spot4toCentreDist','Spot4toCentreDistNormalised',...
    'Spot5X','Spot5Y','Spot5Intensity','Spot5toCentreDist','Spot5toCentreDistNormalised',...
    'Spot6X','Spot6Y','Spot6Intensity','Spot6toCentreDist','Spot6toCentreDistNormalised',...
    'Spot7X','Spot7Y','Spot7Intensity','Spot7toCentreDist','Spot7toCentreDistNormalised',...
    'Spot8X','Spot8Y','Spot8Intensity','Spot8toCentreDist','Spot8toCentreDistNormalised',...
    'Spot9X','Spot9Y','Spot9Intensity','Spot9toCentreDist','Spot9toCentreDistNormalised',...
    'Spot10X','Spot10Y','Spot10Intensity','Spot10toCentreDist','Spot10toCentreDistNormalised',...
    'Spot11X','Spot11Y','Spot11Intensity','Spot11toCentreDist','Spot11toCentreDistNormalised',...
    'Spot12X','Spot12Y','Spot12Intensity','Spot12toCentreDist','Spot12toCentreDistNormalised',...
    'Spot13X','Spot13Y','Spot13Intensity','Spot13toCentreDist','Spot13toCentreDistNormalised',...
    'Spot14X','Spot14Y','Spot14Intensity','Spot14toCentreDist','Spot14toCentreDistNormalised',...
    'Spot15X','Spot15Y','Spot15Intensity','Spot15toCentreDist','Spot15toCentreDistNormalised'};

% save output
writetable(outputTable,[filename,'-Data.txt'])

% plot some graphs and then save these

% plot a distribution of spots

figure;histogram(numSpots')
xlabel('Number of spots in cell');ylabel('Number of cells');
saveas(gcf,strcat(filename,'-barchart-numSpots.png'))

%close(gcf)

% plot the distribution of distances from the centre

for i=1:numel(distances)
    if distances(i) > 0
        nonZeroDists(count) = distances(i);
        count = count + 1;
    end 
end

figure;histogram(nonZeroDists);
xlabel('Distance to centre (pixels)');ylabel('Number of cells');
saveas(gcf,strcat(filename,'-barchart-distancesToCentre.png'))
%close(gcf)

% and normalised distances
nonZeroNormDists = zeros(nnz(normDistances),1);
count = 1;
for i=1:numel(normDistances)
    if normDistances(i) > 0
        nonZeroNormDists(count) = normDistances(i);
        count = count + 1;
    end 
end

figure;histogram(nonZeroNormDists)
xlabel('Normalised distance to centre (pixels)');ylabel('Number of cells');
saveas(gcf,strcat(filename,'-barchart-normDistancesToCentre.png'))

% save a histogram of the number of cells of different radii

figure;histogram(cellRadii,8);
xlabel('Radius (pixels)');ylabel('Number of cells');
saveas(gcf,strcat(filename,'-histogram-radii.png'))
%close(gcf)

close all;
