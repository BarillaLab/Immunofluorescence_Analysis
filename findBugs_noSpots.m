% for Azhar to identify the bugs in her image via DAPI and then use this to
% count how many red spots (a protein?) there are in the image

% outputs a table with columns
% cell number | cell centre x | cell centre y | cell radius...
% Pearson correlation coefficient | minimum red intensity in cell | ...
% maximum red intensity in cell | mean red intensity in cell
% ...

function findBugs_noSpots()

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
cellCentresX = zeros(numAnalysedBugs,1);
cellCentresY = zeros(numAnalysedBugs,1);
cellRadii = zeros(numAnalysedBugs,1);
cellNumbers = (1:numAnalysedBugs)';
PCC = zeros(numAnalysedBugs,1);
distances = zeros(numAnalysedBugs,15);
normDistances = zeros(numAnalysedBugs,15);
minIntensities = zeros(numAnalysedBugs,1);
maxIntensities = zeros(numAnalysedBugs,1);
meanIntensities = zeros(numAnalysedBugs,1);


figure;hold on;
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
    
    plot(redPixelVals,bluePixelVals,'o');
    
    PCC(i) = corr(double(redPixelVals(:)),double(bluePixelVals(:))) ;
    
    % compute the min, max, mean intensity in each red cell
    
    disp(max(redPixelVals(:)));
    disp(max(double(redPixelVals(:))));
    
    max(redCell(:))
    
    minIntensities(i) = min(double(redPixelVals(:)));
    maxIntensities(i) = max(double(redPixelVals(:)));
    meanIntensities(i) = mean(double(redPixelVals(:)));
        
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
    
    % output
    

end


outputTable = table(cellNumbers,cellCentresX,cellCentresY,cellRadii,PCC,...
    minIntensities,maxIntensities,meanIntensities);

% change headings

outputTable.Properties.VariableNames = {'CellNumber','CellCentreX','CellCentreY','CellRadius'...
    'PCC','MinIntensity','MaxIntensity','MeanIntensity'};

% save output
writetable(outputTable,[filename,'-Data.txt'])

% plot some graphs and then save these


% save a histogram of the number of cells of different radii

figure;histogram(cellRadii,8);
xlabel('Radius (pixels)');ylabel('Number of cells');
saveas(gcf,strcat(filename,'-histogram-radii.png'))
%close(gcf)

%close all;
