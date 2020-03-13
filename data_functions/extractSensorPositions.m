function sensorInfo = extractSensorPositions(cfg)
% Extracts the orientations (rad and tan, sensor space)and the position of
% sensors described by STL files.
% 
% This code is quite specific to the STL files for the Gen 2 OPMs, but
% could be adapted for others.
% 
% Inputs. All are required at this stage, but can input empty. 
% cfg.folder        = 'string', folder containing only the STL for sensors.
% cfg.slotori       = 'rad' or 'tan', the sensor-space radial to the head.
% cfg.plot          = bool, do you want a plot of the initial result?
% cfg.correct       = bool, do you want to manually flip the tan ori?
% cfg.scalp         = STL file of scalp in the same coordinate system.
% cfg.outputfile    = 'string', filename with directory for output.
% 
% Output is a table of sensor info. 'rad' and 'tan' are given their own
% channels, rather than combining as Ox_rad, or similar. This should make
% it more compatible with non-dual-axis pipelines. The output is given by
% the function, but also saved to the specied file. Columns as follows:
% 
% 1) filename, STL filename in the order given by dir(). Check it is right.
% 2) slot, the slot number asigned based on the dir() order. 
% 3:5) Px; Py; Pz, position info. Will use whatever units were inputted. 
% 6:8) Ox; Oy; Oz, orientation info. Normalised to the unit inputted. 
% 
% Editted from Stephanie Mellor's (09/02/2020) code by Nicholas Alexander 
% (12/03/2020)

% Input file path and find all STL files. 
STLdir              = strcat(cfg.folder);
sensorListing       = dir([STLdir '\*.stl']);

% We need to know the centre point of all sensors as a rough indication of
% the orientation towards the scalp. Read in all files and take the mean.
% Not the most efficient method!
meanPos         = zeros(size(sensorListing, 1), 3);
for sensorIdx = 1:size(sensorListing, 1)
    sensorFilename          = fullfile(sensorListing(sensorIdx).folder,sensorListing(sensorIdx).name);
    [~, verts]              = stlread(sensorFilename);
    meanPos(sensorIdx,1:3)  = mean(verts,1);
end
headCentreVert  = mean(meanPos,1);

clear meanPos sensorIdx verts sensorFilename STLdir


% Preallocate the main outputs. 
centrePoint         = zeros(size(sensorListing, 1), 3);
radOri              = zeros(size(sensorListing, 1), 3);
tanOri              = zeros(size(sensorListing, 1), 3);
sensorFaces         = cell(size(sensorListing, 1), 3);
sensorVerts         = cell(size(sensorListing, 1), 3);
filename_rad        = cell(size(sensorListing, 1), 3);
filename_tan        = cell(size(sensorListing, 1), 3);

for sensorIdx = 1:size(sensorListing, 1)
    % Specifiy and read in sensor STLs one at a time. 
    sensorFilename  = fullfile(sensorListing(sensorIdx).folder,sensorListing(sensorIdx).name);
    [faces, verts]              = stlread(sensorFilename);
    
    filename_rad{sensorIdx}     = strcat(sensorFilename,'_rad');
    filename_tan{sensorIdx}     = strcat(sensorFilename,'_tan');
    
    % Reduce the vertices by removing duplicate positions. This will help
    % later on. For example, it avoids there being duplicate corners which
    % would prevent using rules like only having 8 corners. 
    [redVerts,~,~]      = unique(verts,'rows');
    
    % Face references need to be updated to match them.
    redFaces            = zeros(size(faces));
    for redVertIdx = 1:length(redVerts(:,1))
        vertRepeatIdx           = find(ismember(verts,redVerts(redVertIdx,:),'rows'));
        faceRepeatIdx           = find(ismember(faces,vertRepeatIdx));
        redFaces(faceRepeatIdx) = redVertIdx; %#ok<FNDSB>
    end
    
    % We can find the corners by finding the points furthest from the mean
    % of all vertices. Note, this may not be the centre.
    meanVert        = mean(redVerts);
    meanVertDist    = zeros(length(redVerts(:,1)),1);
    for pointIdx = 1:length(redVerts(:,1))
        meanVertDist(pointIdx)  = pdist([meanVert;redVerts(pointIdx,1:3)]);
    end
    
    % Just the 8 greatest distances.
    [~,distRank]        = sort(meanVertDist);
    cornerIdx           = distRank(end-7:end);
    cornerVerts         = redVerts(cornerIdx,1:3);
    
    % Up to here, everything works. But need a different way of finding the
    % groups of sensors (top/bottom).
    
    % Find the distance between all the corner combinations.
    combinations        = nchoosek(1:length(cornerVerts(:,1)),2);
    cornerDist          = zeros(length(combinations(:,1)),1);
    for combIdx = 1:length(combinations(:,1))
        cornerDist(combIdx,1)   = pdist([cornerVerts(combinations(combIdx,1),:);cornerVerts(combinations(combIdx,2),:)]);
    end
    
    % Round it out and sort it out.
    cornerDist          = round(cornerDist,2);
    uniqueCornerDists   = sort(unique(cornerDist));
    
    % Remove the hypotenusesese.
    combinations    = nchoosek(1:length(uniqueCornerDists),2);
    hypoIdx       = zeros(length(uniqueCornerDists),length(combinations(:,1)));
    for combIdx = 1:length(combinations(:,1))
        RMS             = sqrt(uniqueCornerDists(combinations(combIdx,1))^2 + uniqueCornerDists(combinations(combIdx,2))^2);
        RMS             = round(RMS,2);
        [~,hypoIdx(:,combIdx)]    = ismember(uniqueCornerDists,RMS);
    end
    hypoIdx     = (sum(hypoIdx,2)>0);
    
    % The remaining distances are along the outside edges.
    edgeDists   = uniqueCornerDists(~hypoIdx);
    
    % If the slot rad and head rad are aligned than it is short edge down
    % so remove the longest edge. Otherwise it is long edge, so remove the
    % shortest edge.
    if strcmp(cfg.slotori,'rad')
        edgeDists(edgeDists == max(edgeDists)) = [];
    elseif strcmp(cfg.slotori,'tan')
        edgeDists(edgeDists == min(edgeDists)) = [];
    else
        error('Specify either cfg.slotori = rad or tan')
    end
    
    % Take a vertex and find the connected vertices. 
    firstGroupCornerVerts   = zeros(length(cornerVerts(:,1))/2,length(cornerVerts(1,:)));
    remCornerVerts          = cornerVerts;
    nextDist                = 1;
    for cornerIdx = 1:length(firstGroupCornerVerts(:,1))
        if cornerIdx == 1
            firstGroupCornerVerts(cornerIdx,:)  = remCornerVerts(cornerIdx,:);
            remCornerVerts(1,:)                 = [];
        else
            dist        = zeros(length(remCornerVerts(:,1)),1);
            for remIdx = 1:length(remCornerVerts(:,1))
                dist(remIdx,1)  = pdist([firstGroupCornerVerts(cornerIdx-1,:);remCornerVerts(remIdx,:)]);
            end
            dist    = round(dist,2);
            % Long or short edge this time?
            if nextDist == 1
                nextDist = 2;
            elseif nextDist == 2
                nextDist = 1;
            end

            nextCornerIdx                       = (dist == edgeDists(nextDist));
            firstGroupCornerVerts(cornerIdx,:)  = remCornerVerts(nextCornerIdx,:);
            remCornerVerts(nextCornerIdx,:)         = [];
        end
    end
    
    % The remaining corner verts are the second group.
    secondGroupCornerVerts  = remCornerVerts;
    
    % Take the mean position of each grouping.
    firstGroupCentreVert    = mean(firstGroupCornerVerts,1);
    secondGroupCentreVert   = mean(secondGroupCornerVerts,1);
    
    % Which is closer to the centre of the head?
    firstGroupCentreDist    = pdist([firstGroupCentreVert;headCentreVert]);
    secondGroupCentreDist   = pdist([secondGroupCentreVert;headCentreVert]);
    
    if (firstGroupCentreDist < secondGroupCentreDist)
        bottomCornerVerts   = firstGroupCornerVerts;
        topCornerVerts      = secondGroupCornerVerts;
    elseif (firstGroupCentreDist > secondGroupCentreDist)
        bottomCornerVerts   = secondGroupCornerVerts;
        topCornerVerts      = firstGroupCornerVerts;
    end   

    % Find the distances between the top corner points.
    combinations = nchoosek(1:4, 2);
    topCornerDistance    = zeros(length(combinations),1);
    bottomCornerDistance = zeros(length(combinations),1);
    for combIdx = 1:length(combinations)
        topCornerDistance(combIdx)      = pdist(vertcat(topCornerVerts(combinations(combIdx,1),:),topCornerVerts(combinations(combIdx,2),:)));
        bottomCornerDistance(combIdx)   = pdist(vertcat(bottomCornerVerts(combinations(combIdx,1),:),bottomCornerVerts(combinations(combIdx,2),:)));
    end
    
    % Longest distance goes across the middle.
    [~,topMiddleComb]           = max(topCornerDistance);
    [~,bottomMiddleComb]        = max(bottomCornerDistance);

    % Shortest distance goes along the ends.
    [~, bottomEndComb]          = min(bottomCornerDistance);
    
    % Find the centre point. Average position of the furthest apart pairs. 
    topMiddlePoint              = mean(vertcat(topCornerVerts(combinations(topMiddleComb,1),:),topCornerVerts(combinations(topMiddleComb,2),:)));
    bottomMiddlePoint           = mean(vertcat(bottomCornerVerts(combinations(bottomMiddleComb,1),:),bottomCornerVerts(combinations(bottomMiddleComb,2),:)));
    centrePoint(sensorIdx,1:3)  = mean(vertcat(topMiddlePoint,bottomMiddlePoint));   
    
    bottomEndPoint              = mean(vertcat(bottomCornerVerts(combinations(bottomEndComb,1),:),bottomCornerVerts(combinations(bottomEndComb,2),:)));
    
    % Very confusing, but if the sensors are in 'tan' orientation slots,
    % then that means the tan of the sensor is aligned with the rad of the
    % head. Need to avvound for that. 
    if strcmp(cfg.slotori,'rad')
        radOri(sensorIdx,1:3)       = bottomMiddlePoint - centrePoint(sensorIdx,1:3);
        tanOri(sensorIdx,1:3)       = bottomMiddlePoint - bottomEndPoint(1:3);
    elseif strcmp(cfg.slotori,'tan')
        radOri(sensorIdx,1:3)       = bottomMiddlePoint - centrePoint(sensorIdx,1:3);
        tanOri(sensorIdx,1:3)       = bottomMiddlePoint - bottomEndPoint(1:3);
    end
    
    % Normalise the orientations
    L                           = sqrt(radOri(sensorIdx,1)^2 + radOri(sensorIdx,2)^2 + radOri(sensorIdx,3)^2);
    radOri(sensorIdx,1:3)       = radOri(sensorIdx,1:3)/L;
    L                           = sqrt(tanOri(sensorIdx,1)^2 + tanOri(sensorIdx,2)^2 + tanOri(sensorIdx,3)^2);
    tanOri(sensorIdx,1:3)       = tanOri(sensorIdx,1:3)/L;
    
    % Make a matlab mesh structure from the reduced faces and verts.
    sensorFaces{sensorIdx}(:,1:3)   = redFaces;
    sensorVerts{sensorIdx}(:,1:3)   = redVerts;
end

% Plot the result
if cfg.plot
    figure(1);
    hold on;
    grid on;
    quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), radOri(:,1), radOri(:,2), radOri(:,3));
    quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), tanOri(:,1), tanOri(:,2), tanOri(:,3));
    scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
    text(centrePoint(:,1)-1.5*radOri(:,1),centrePoint(:,2)-1.5*radOri(:,2),centrePoint(:,3)-1.5*radOri(:,3), cellstr(num2str((1:size(centrePoint,1))')), 'color', 'g');
    daspect([1 1 1])
    for sensorIdx = 1:size(sensorListing, 1)
        patch('Faces',sensorFaces{sensorIdx},'Vertices',sensorVerts{sensorIdx},'FaceAlpha',.3,'EdgeAlpha',0);
    end
    
    if ~isempty(cfg.scalp)
        [scalpFaces, scalpVerts]              = stlread(cfg.scalp);
        patch('Faces',scalpFaces,'Vertices',scalpVerts,'FaceAlpha',.1,'EdgeAlpha',0);
    end
    hold off
    
    % Ask if they are done with it.
    input('Press any key to continue (closes figure)')
    close all % Change this later. 
end

% For determining whether to flip the orientation or not plot one sensor
% onto the scalp at a time.
if cfg.correct
    repeat = 1;
    while repeat
        if ~isempty(cfg.scalp)
                [scalpFaces, scalpVerts]              = stlread(cfg.scalp);
        end
        for sensorIdx = 1:size(sensorListing, 1)
            hold on
            if ~isempty(cfg.scalp)
                patch('Faces',scalpFaces,'Vertices',scalpVerts,'FaceAlpha',1,'EdgeAlpha',0,'FaceColor',[.85 .72 .6]);
            end
            patch('Faces',sensorFaces{sensorIdx},'Vertices',sensorVerts{sensorIdx},'FaceAlpha',1,'EdgeAlpha',0,'FaceColor',[0 0 0])
            quiver3(centrePoint(sensorIdx,1), centrePoint(sensorIdx,2), centrePoint(sensorIdx,3), tanOri(sensorIdx,1).*50, tanOri(sensorIdx,2).*50, tanOri(sensorIdx,3).*50,'color',[0 0 1]);
            scatter3(centrePoint(sensorIdx,1), centrePoint(sensorIdx,2), centrePoint(sensorIdx,3))
            view(-radOri(sensorIdx,:));
            
            flip        = input('Press 1 to flip or 0 to skip');
            if flip == 1
                tanOri(sensorIdx,1:3)   = -tanOri(sensorIdx,1:3);
            end
            if sensorIdx ~= size(sensorListing, 1)
                close all
            end
        end
        repeat    = input('Would you like to go through them again? 1 for yes, 0 for no');
        close all
    end
end

% Put the main info into one output.
filename    = [filename_rad;filename_tan];
slot        = [1:size(sensorListing, 1);1:size(sensorListing, 1);];
Px          = [centrePoint(:,1);centrePoint(:,1)];
Py          = [centrePoint(:,2);centrePoint(:,2)];
Pz          = [centrePoint(:,3);centrePoint(:,3)];
Ox          = [radOri(:,1);tanOri(:,1)];
Oy          = [radOri(:,2);tanOri(:,2)];
Oz          = [radOri(:,3);tanOri(:,3)];
outputTable = table(filename,slot,Px,Py,Pz,Ox,Oy,Oz);
if ~isempty(cfg.outputfile)
    
    writetable(outputTable,strcat(cfg.outputfile,'.tsv'),'Delimiter','tab','QuoteStrings',false,'FileType', 'text');
end

sensorInfo      = outputTable;