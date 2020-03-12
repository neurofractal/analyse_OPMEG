function grad = extractSensorPositions(folder,slotOri)
% Extracts the orientations (rad and tan, sensor space)and the position of
% sensors described by STL files.
%
% Input a folder containing only STL files for individual sensors. They
% should all be in the same coordinate system. Describe the slot
% orientation as 'rad' or 'tan', indicating which orientation at the sensor
% level is pointing towards the scalp. 

% Input file path and find all STL files. 
STLdir              = strcat(folder);
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

clear meanPos sensorIdx verts sensorFilename STLdir folder


% Preallocate the main outputs. 
centrePoint         = zeros(size(sensorListing, 1), 3);
radOri              = zeros(size(sensorListing, 1), 3);
tanOri              = zeros(size(sensorListing, 1), 3);
sensorFaces         = cell(size(sensorListing, 1), 3);
sensorVerts         = cell(size(sensorListing, 1), 3);

for sensorIdx = 1:size(sensorListing, 1)
    % Specifiy and read in sensor STLs one at a time. 
    sensorFilename      = fullfile(sensorListing(sensorIdx).folder,sensorListing(sensorIdx).name);
    [faces, verts]      = stlread(sensorFilename);
    
    % Reduce the vertices by removing duplicate positions. This will help
    % later on. For example, it avoids there being duplicate corners which
    % would prevent using rules like only having 8 corners. 
    [redVerts,~,~]      = unique(verts,'rows');
    
    % Face references need to be updated to match them.
    redFaces            = zeros(size(faces));
    for redVertIdx = 1:length(redVerts(:,1))
        vertRepeatIdx           = find(ismember(verts,redVerts(redVertIdx,:),'rows'));
        faceRepeatIdx           = find(ismember(faces,vertRepeatIdx));
        redFaces(faceRepeatIdx) = redVertIdx;
    end
    
    % We can find the corners by finding the points furthest from the mean
    % of all vertices. Note, this may not be the centre.
    meanVert        = mean(redVerts);
    meanVertDist    = zeros(length(redVerts(:,1)),1);
    for pointIdx = 1:length(redVerts(:,1))
        meanVertDist(pointIdx)  = pdist([meanVert;redVerts(pointIdx,1:3)]);
    end
    
    % Just the 8 greatest distances.
    [~,distRank]         = sort(meanVertDist);
    cornerIdx           = distRank(end-7:end);
    cornerVerts         = redVerts(cornerIdx,1:3);
    
    % Find the distance between each sensor and the centre of the head.
    cornerCentreDist    = zeros(length(cornerVerts(:,1)),1);
    for cornerIdx = 1:length(cornerVerts(:,1))
        cornerCentreDist(cornerIdx,1) = pdist([cornerVerts(cornerIdx,1:3);headCentreVert]);
    end
    
    % The closest four will be the bottom corners.
    bottomCornerIdx     = (cornerCentreDist < median(cornerCentreDist));
    bottomCornerVerts   = cornerVerts(bottomCornerIdx,1:3);
    topCornerVerts      = cornerVerts(~bottomCornerIdx,1:3);

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
    
    % Add orientation (tan)
    radOri(sensorIdx,1:3)       = bottomMiddlePoint - centrePoint(sensorIdx,1:3);
        
    % Add orientation (rad)
    tanOri(sensorIdx,1:3)       = bottomMiddlePoint - bottomEndPoint(1:3);
    
    % Normalise the orientations
    L                           = sqrt(radOri(sensorIdx,1)^2 + radOri(sensorIdx,2)^2 + radOri(sensorIdx,3)^2);
    radOri(sensorIdx,1:3)       = radOri(sensorIdx,1:3)/L;
    L                           = sqrt(tanOri(sensorIdx,1)^2 + tanOri(sensorIdx,2)^2 + tanOri(sensorIdx,3)^2);
    tanOri(sensorIdx,1:3)       = tanOri(sensorIdx,1:3)/L;
    
    % Make a matlab mesh structure from the reduced faces and verts.
    sensorFaces{sensorIdx}(:,1:3)   = redFaces;
    sensorVerts{sensorIdx}(:,1:3)   = redVerts;
end

%% Show
figure;
hold on;
grid on;
quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), radOri(:,1), radOri(:,2), radOri(:,3));
quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), tanOri(:,1), tanOri(:,2), tanOri(:,3));

scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
text(centrePoint(:,1)-1.5*radOri(:,1),centrePoint(:,2)-1.5*radOri(:,2),centrePoint(:,3)-1.5*radOri(:,3), cellstr(num2str((1:size(centrePoint,1))')), 'color', 'g');
daspect([1 1 1])

% plot the STLs as well
for sensorIdx = 1:size(sensorListing, 1)
    patch('Faces',sensorFaces{sensorIdx},'Vertices',sensorVerts{sensorIdx},'FaceAlpha',.1,'EdgeAlpha',0);

end
hold off

% 
% % For determining whether to flip the orientation or not plot one sensor
% % onto the scalp at a time.
% STLdir2 = strcat(workingDir,'\Data\STL_Positions\Hingecast\Gareth Head.stl');
% sensorListing2 = dir(STLdir2);
% [F3, V3] = stlread(fullfile(sensorListing2.folder, sensorListing2.name));
% 
% 
% % plot the STLs as well
% for sensorIdx = 1:size(sensorListing, 1)
%     hold on
%     patch('Faces',Mesh.Faces{sensorIdx},'Vertices',Mesh.Points{sensorIdx},'FaceAlpha',.1,'EdgeAlpha',0)
%     quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), tanOri(:,1), tanOri(:,2), tanOri(:,3));
%     scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
%     flip = input('flip?');
%     if flip == 1
%         decision(sensorIdx) = true;
%         tanOri(sensorIdx,1:3)    = -tanOri(sensorIdx,1:3);
%     else
%         decision(sensorIdx) = false;
%     end
%     close all
% end
% 
% % Alternatively, specify the numbers to be flipped.
% flipSensors = input('Form [1 56]: ');
% 
% for sensorIdx = 1:length(flipSensors)
%     tanOri(flipSensors(sensorIdx),1:3)    = -tanOri(flipSensors(sensorIdx),1:3);
% end
% 
% 
% 
% %% Show
% figure;
% hold on;
% grid on;
% quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), normZrad(:,1), normZrad(:,2), normZrad(:,3));
% quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), normZtan(:,1), normZtan(:,2), normZtan(:,3));
% 
% scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
% text(centrePoint(:,1)-1.5*normZrad(:,1),centrePoint(:,2)-1.5*normZrad(:,2),centrePoint(:,3)-1.5*normZrad(:,3), cellstr(num2str((1:size(centrePoint,1))')), 'color', 'g');
% daspect([1 1 1])
% 
% % plot the STLs as well
% for sensorIdx = 1:size(sensorListing, 1)
%     patch('Faces',Mesh.Faces{sensorIdx},'Vertices',Mesh.Points{sensorIdx},'FaceAlpha',.1,'EdgeAlpha',0)
% end
% 
% patch('Faces',F3,'Vertices',V3,'FaceAlpha',.1,'EdgeAlpha',0)
% 
% hold off
% 
% % Make the grad structure.
% for sensorIdx = 1:size(sensorListing, 1)
%     label{sensorIdx}    = strcat(num2str(sensorIdx),'_rad');
%     label{sensorIdx+size(sensorListing, 1)}     = strcat(num2str(sensorIdx),'_tan');
% end
% 
% grad            = [];
% grad.chanpos    = [centrePoint; centrePoint];
% grad.chanori    = [normZrad; normZtan];
% 
% 
% grad.label      = label';
% 
% 
% %% Next step is to combine this with a file to get the actual position of the cell. 
% % Provide some info on which way the sensor is put into the slot. 
% 
