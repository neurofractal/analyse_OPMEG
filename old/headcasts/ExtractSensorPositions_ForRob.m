
% TODO:
% - fix axis problem - if sensor is perfectly aligned with one axis, can't
% find the corner points the way I've currently done it
% - find a way to set front of head for tangential sensors. For the radial
% ones, choose the sensor end either from the points on Mark's STL file or
% the center of the headcast
%
% NA Edit of Steph's code

%% Setup Directories.
% pcName          = getenv('COMPUTERNAME');
% if strcmp(pcName,'PRINCE')
%     matlabDir       = 'D:\MATLAB';
% else
%     disp('Will need to set matlabDir manually. Probably the rest too');
% end
% 
% toolboxesDir     = strcat(matlabDir,'\Toolboxes');
% fieldtripDir    = strcat(toolboxesDir,'\fieldtrip-20190419');
% workingDir      = strcat(matlabDir,'\Analysis\OPM_Analysis');
% functionsDir    = genpath(strcat(workingDir,'\Functions'));
% templatesDir    = strcat(workingDir,'\Templates');
% 
% addpath(fieldtripDir)
% ft_defaults;
% cd(workingDir)
% addpath(functionsDir,templatesDir,workingDir)
% 
% clear templatesDir fieldtripDir matlabDir functionsDir atlasDir pcName

%% config
% N.B. if sensors are lying flat, set sensorPoints to false. The sensitive
% region will be set as the center of the OPM, rather than 6 mm
sensorsRad = false;

% Input file path
STLdir = 'D:\data\gareth_scanner_cast\indiv_stl';

% Check if sensor stls have points at sensor end
sensorPoints = true;

%% main
sensorListing = dir([STLdir '\*.stl']);
centrePoint = zeros(size(sensorListing, 1), 3);
% X = zeros(size(sensorListing, 1), 3);
% Y = zeros(size(sensorListing, 1), 3);
Z = zeros(size(sensorListing, 1), 3);
% 
% if ~sensorPoints
%     faceCentersVec = zeros(2, 3, size(sensorListing,1));
% end

for i = 1:size(sensorListing, 1)
    if i == 34
        disp('');
    end
    
    [F, V] = stlread(fullfile(sensorListing(i).folder, sensorListing(i).name));
    
    % Remove duplicate positions.
    [V2,~,~] = unique(V,'rows');
    F2              = zeros(size(F));
    for uniqueIdx = 1:length(V2(:,1))
        pointRepeatIdx      = find(ismember(V,V2(uniqueIdx,:),'rows'));
        pointRepeatIdx(1);
        faceRepeatIdx       = find(ismember(F,pointRepeatIdx));
        F2(faceRepeatIdx)   = uniqueIdx;
    end
    
    coneIdx    = mode(F2,'all');
    conePoint   = V2(coneIdx,:);
    
    % Find index of 1st corner point
    [~, I{1}] = min(V(:,1), [], 1);
    
    % Find index of 2nd corner point
    [~, I{2}] = max(V(:,1), [], 1);
    
    % Find index 3rd corner
    [~, I{3}] = min(V(:,2), [], 1);
    
    % Find index 4th corner
    [~, I{4}] = max(V(:,2), [], 1);
    
    % Find 5th corner
    [~, I{5}] = min(V(:,3), [], 1);
    
    % Find 6th corner
    [~, I{6}] = max(V(:,3), [], 1);
    
    for j = 1:6
        if size(I{j}) > 1
            disp('Error - Sensor is aligned with axis');
        end
    end
    
    % Find sides between corners
    combinations = nchoosek(1:6, 2);
    s = zeros(size(combinations, 1), 3);
    for j = 1:size(combinations, 1)
        s(j, :) = V(I{combinations(j,1)},:) - V(I{combinations(j,2)},:);
    end
   
    cornerPoints     = zeros(6,3);
    for cornerIdx = 1:6
        cornerPoints(cornerIdx,1:3)  = V(I{cornerIdx},1:3);
    end
    
    
    
    % Distance from each corner to the cone. 
    cornerConeDist = zeros(1,6);
    for cornerIdx = 1:6
        cornerConeDist(cornerIdx)  = pdist(vertcat(V2(coneIdx,:),cornerPoints(cornerIdx,:)));
    end
    
    % The 3 greater cornerCondDist are the top and the three lesser are the
    % bottom. Identify them.
    topCornerIdx        = (cornerConeDist > median(cornerConeDist));
    bottomCornerIdx     = (cornerConeDist < median(cornerConeDist));
    
    topCornerPoints     = cornerPoints(topCornerIdx,:);
    bottomCornerPoints  = cornerPoints(bottomCornerIdx,:);
    
    % Find the distances between the top corner points.
    combinations = nchoosek(1:3, 2);
    for combIdx = 1:length(combinations)
        topCornerDistance(combIdx)      = pdist(vertcat(topCornerPoints(combinations(combIdx,1),:),topCornerPoints(combinations(combIdx,2),:)));
        bottomCornerDistance(combIdx)   = pdist(vertcat(bottomCornerPoints(combinations(combIdx,1),:),bottomCornerPoints(combinations(combIdx,2),:)));
    end
    
    % Longest distance goes across the middle.
    [~,topMiddleComb] = max(topCornerDistance);
    [~,bottomMiddleComb] = max(bottomCornerDistance);

    % Starting to realise top and bottom are the same. At least for this
    % one. Oh well. Also, bottom might not necessarily be bottom. Use cone
    % later to check. 
    
    % Average position of the furthest apart pairs. 
    topMiddlePoint          = mean(vertcat(topCornerPoints(combinations(topMiddleComb,1),:),topCornerPoints(combinations(topMiddleComb,2),:)));
    bottomMiddlePoint       = mean(vertcat(bottomCornerPoints(combinations(bottomMiddleComb,1),:),bottomCornerPoints(combinations(bottomMiddleComb,2),:)));
    centrePoint(i,1:3)      = mean(vertcat(topMiddlePoint,bottomMiddlePoint));   
end
    
%% Show
figure;
hold on;
grid on;
quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), Z(:,1), Z(:,2), Z(:,3));
scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
text(centrePoint(:,1)-1.5*Z(:,1),centrePoint(:,2)-1.5*Z(:,2),centrePoint(:,3)-1.5*Z(:,3), cellstr(num2str((1:size(centrePoint,1))')), 'color', 'g');
daspect([1 1 1])

% plot the STLs as well
for i = 1:size(sensorListing, 1)
    patch('Faces',Mesh.Faces{i},'Vertices',Mesh.Points{i},'FaceAlpha',.1,'EdgeAlpha',.25)
end
hold off