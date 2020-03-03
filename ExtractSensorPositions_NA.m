% Script to save all sensor locations and orientations for headcasts
%
% NA Edit of Steph's code

%% Setup Directories.
pcName          = getenv('COMPUTERNAME');
if strcmp(pcName,'PRINCE')
    matlabDir       = 'D:\MATLAB';
else
    disp('Will need to set matlabDir manually. Probably the rest too');
end

toolboxesDir     = strcat(matlabDir,'\Toolboxes');
fieldtripDir    = strcat(toolboxesDir,'\fieldtrip-20190419');
workingDir      = strcat(matlabDir,'\Analysis\OPM_Analysis');
functionsDir    = genpath(strcat(workingDir,'\Functions'));
templatesDir    = strcat(workingDir,'\Templates');

addpath(fieldtripDir)
ft_defaults;
cd(workingDir)
addpath(functionsDir,templatesDir,workingDir)

clear templatesDir fieldtripDir matlabDir functionsDir atlasDir pcName

%% config
% N.B. if sensors are lying flat, set sensorPoints to false. The sensitive
% region will be set as the center of the OPM, rather than 6 mm
sensorsRad = false;

% Input file path
STLdir = strcat(workingDir,'\Data\STL_Positions\Hingecast\SensorSTLs_RS');

% Check if sensor stls have points at sensor end
sensorPoints = true;

%% main
sensorListing = dir([STLdir '\*.stl']);
centrePoint = zeros(size(sensorListing, 1), 3);
Zrad    = zeros(size(sensorListing, 1), 3);
Ztan    = zeros(size(sensorListing, 1), 3);

for i = 1:size(sensorListing, 1)
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
    [~, I{1}] = min(V2(:,1), [], 1);
    
    % Find index of 2nd corner point
    [~, I{2}] = max(V2(:,1), [], 1);
    
    % Find index 3rd corner
    [~, I{3}] = min(V2(:,2), [], 1);
    
    % Find index 4th corner
    [~, I{4}] = max(V2(:,2), [], 1);
    
    % Find 5th corner
    [~, I{5}] = min(V2(:,3), [], 1);
    
    % Find 6th corner
    [~, I{6}] = max(V2(:,3), [], 1);
    
    % Are all the corners unique?
    fixAlignedSensor = false;
    if length(unique(cell2mat(I))) < 6
        disp('Error - Sensor is aligned with axis');
        fixAlignedSensor = true;
    end
    
    if fixAlignedSensor
        % How else can I find the corners?
        % Furthest point from cone to get top.
        for pointIdx = 1:length(V2(:,1))
            pointConeDist(pointIdx) = pdist([conePoint;V2(pointIdx,1:3)]);
        end
        pointConeDist = round(pointConeDist,3);
        topCornerPointsIdx = find(pointConeDist == max(pointConeDist));
        topCornerPoints     = V2(topCornerPointsIdx,1:3);
        % Find the furthest away points from each of those. Should be 4. To
        % work with the rest of the code just remove one of them.
        topCornerPoints = topCornerPoints(1:3,1:3);
        for topCornerPointsIdx = 1:length(topCornerPoints(:,1))
            for pointIdx = 1:length(V2(:,1))
                pointCornerDist(pointIdx) = pdist([topCornerPoints(topCornerPointsIdx,1:3);V2(pointIdx,1:3)]);
            end
            pointCornerDist = round(pointCornerDist,3);
            % there should just be one point that is far away. The opposite
            % corner.
            bottomCornerPointsIdx(topCornerPointsIdx) = find(pointCornerDist == max(pointCornerDist));
        end
        bottomCornerPoints     = V2(bottomCornerPointsIdx,1:3);
    else


        % Find sides between corners
        combinations = nchoosek(1:6, 2);
        s = zeros(size(combinations, 1), 3);
        for j = 1:size(combinations, 1)
            s(j, :) = V2(I{combinations(j,1)},:) - V2(I{combinations(j,2)},:);
        end

        cornerPoints     = zeros(6,3);
        for cornerIdx = 1:6
            cornerPoints(cornerIdx,1:3)  = V2(I{cornerIdx},1:3);
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
    end

    % Find the distances between the top corner points.
    combinations = nchoosek(1:3, 2);
    for combIdx = 1:length(combinations)
        topCornerDistance(combIdx)      = pdist(vertcat(topCornerPoints(combinations(combIdx,1),:),topCornerPoints(combinations(combIdx,2),:)));
        bottomCornerDistance(combIdx)   = pdist(vertcat(bottomCornerPoints(combinations(combIdx,1),:),bottomCornerPoints(combinations(combIdx,2),:)));
    end
    
    % Longest distance goes across the middle.
    [~,topMiddleComb]       = max(topCornerDistance);
    [~,bottomMiddleComb]    = max(bottomCornerDistance);

    % Shortest distance goes along the ends.
    [~, bottomEndComb]      = min(bottomCornerDistance);
    
    % Average position of the furthest apart pairs. 
    topMiddlePoint          = mean(vertcat(topCornerPoints(combinations(topMiddleComb,1),:),topCornerPoints(combinations(topMiddleComb,2),:)));
    bottomMiddlePoint       = mean(vertcat(bottomCornerPoints(combinations(bottomMiddleComb,1),:),bottomCornerPoints(combinations(bottomMiddleComb,2),:)));
    centrePoint(i,1:3)      = mean(vertcat(topMiddlePoint,bottomMiddlePoint));   
    
    bottomEndPoint          = mean(vertcat(bottomCornerPoints(combinations(bottomEndComb,1),:),bottomCornerPoints(combinations(bottomEndComb,2),:)));

    
    % Add direction (tan)
    Zrad(i,1:3)                = conePoint - centrePoint(i,1:3);
    Mesh.Faces{i}(:,1:3)       = F2;
    Mesh.Points{i}(:,1:3)      = V2;
        
    % Add direction (rad)
    Ztan(i,1:3)                 = conePoint - bottomEndPoint(1:3);
    
end

%% Show
figure;
hold on;
grid on;
quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), Zrad(:,1), Zrad(:,2), Zrad(:,3));
quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), Ztan(:,1), Ztan(:,2), Ztan(:,3));

scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
text(centrePoint(:,1)-1.5*Zrad(:,1),centrePoint(:,2)-1.5*Zrad(:,2),centrePoint(:,3)-1.5*Zrad(:,3), cellstr(num2str((1:size(centrePoint,1))')), 'color', 'g');
daspect([1 1 1])

% plot the STLs as well
for i = 1:size(sensorListing, 1)
    patch('Faces',Mesh.Faces{i},'Vertices',Mesh.Points{i},'FaceAlpha',.1,'EdgeAlpha',0)
end
hold off

% For determining whether to flip the orientation or not plot one sensor
% onto the scalp at a time.
STLdir2 = strcat(workingDir,'\Data\STL_Positions\Hingecast\Gareth Head.stl');
sensorListing2 = dir(STLdir2);
[F3, V3] = stlread(fullfile(sensorListing2.folder, sensorListing2.name));


% plot the STLs as well
for i = 1:size(sensorListing, 1)
    hold on
    patch('Faces',Mesh.Faces{i},'Vertices',Mesh.Points{i},'FaceAlpha',.1,'EdgeAlpha',0)
    quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), Ztan(:,1), Ztan(:,2), Ztan(:,3));
    scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
    flip = input('flip?');
    if flip == 1
        decision(i) = true;
        Ztan(i,1:3)    = -Ztan(i,1:3);
    else
        decision(i) = false;
    end
    close all
end

% Alternatively, specify the numbers to be flipped.
flipSensors = input('Form [1 56]: ');

for i = 1:length(flipSensors)
    Ztan(flipSensors(i),1:3)    = -Ztan(flipSensors(i),1:3);
end



%% Normalise
for i = 1:length(Zrad(:,1))
    L               = sqrt(Zrad(i,1)^2 + Zrad(i,2)^2 + Zrad(i,3)^2);
    normZrad(i,1:3) = Zrad(i,1:3)/L;
    L               = sqrt(Ztan(i,1)^2 + Ztan(i,2)^2 + Ztan(i,3)^2);
    normZtan(i,1:3) = Ztan(i,1:3)/L;
end

%% Show
figure;
hold on;
grid on;
quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), normZrad(:,1), normZrad(:,2), normZrad(:,3));
quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3), normZtan(:,1), normZtan(:,2), normZtan(:,3));

scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
text(centrePoint(:,1)-1.5*normZrad(:,1),centrePoint(:,2)-1.5*normZrad(:,2),centrePoint(:,3)-1.5*normZrad(:,3), cellstr(num2str((1:size(centrePoint,1))')), 'color', 'g');
daspect([1 1 1])

% plot the STLs as well
for i = 1:size(sensorListing, 1)
    patch('Faces',Mesh.Faces{i},'Vertices',Mesh.Points{i},'FaceAlpha',.1,'EdgeAlpha',0)
end

patch('Faces',F3,'Vertices',V3,'FaceAlpha',.1,'EdgeAlpha',0)

hold off

% Make the grad structure.
for i = 1:size(sensorListing, 1)
    label{i}    = strcat(num2str(i),'_rad');
    label{i+size(sensorListing, 1)}     = strcat(num2str(i),'_tan');
end

grad            = [];
grad.chanpos    = [centrePoint; centrePoint];
grad.chanori    = [normZrad; normZtan];


grad.label      = label';


%% Next step is to combine this with a file to get the actual position of the cell. 
% Provide some info on which way the sensor is put into the slot. 

