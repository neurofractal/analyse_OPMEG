% Script to save all sensor locations and orientations for headcasts

% TODO:
% - fix axis problem - if sensor is perfectly aligned with one axis, can't
% find the corner points the way I've currently done it
% - find a way to set front of head for tangential sensors. For the radial
% ones, choose the sensor end either from the points on Mark's STL file or
% the center of the headcast

%% Housekeeping
clear;
close all;

%% config

% N.B. if sensors are lying flat, set sensorPoints to false. The sensitive
% region will be set as the center of the OPM, rather than 6 mm
sensorsRad = false;

% Input file path
sensorSTLs = 'C:\Users\rseymour\Dropbox\Research\Projects\2020\headcasts\Share\Hingecast\SensorSTLs\';
% Output file path
savePath = 'D:\Github\analyse_OPMEG\headcasts\';

% Check if sensor stls have points at sensor end
sensorPoints = false;

%% main

sensorListing = dir([sensorSTLs '/*.stl']);
pos = zeros(size(sensorListing, 1), 3);
X = zeros(size(sensorListing, 1), 3);
Y = zeros(size(sensorListing, 1), 3);
Z = zeros(size(sensorListing, 1), 3);

if ~sensorPoints
    faceCentersVec = zeros(2, 3, size(sensorListing,1));
end

for i = 1:size(sensorListing, 1)
    if i == 34
        disp('');
    end
    
    [~, V, ~] = stlread(fullfile(sensorListing(i).folder, sensorListing(i).name));
    
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
   
   % Find shortest and next shortest sides
   sideLengths = sqrt(sum(s.^2, 2));
   
   [shortestSideLength, shortestSideIdx] = min(sideLengths);
   sideLengths(abs(sideLengths - shortestSideLength) < 1e-3) = max(sideLengths);
   [nextshortestSideLength, nextshortestSideIdx] = min(sideLengths);
   
   % Find directions of axes (give or take a factor of -1)
   Ytmp = s(shortestSideIdx, :)./shortestSideLength;
   X(i,:) = s(nextshortestSideIdx, :)./nextshortestSideLength;
   Ztmp = cross(X(i,:), Ytmp);
   Ztmp = Ztmp./norm(Ztmp);
   
   % Find center of cube
   tmp = cell2mat(I);
   vertices = V(tmp,:);
   clear tmp
   center = mean(vertices, 1);

   % Find centers of smaller faces
   if sensorPoints
       faceCenters(1,:) = center - Ztmp*24/2;
       faceCenters(2,:) = center + Ztmp*24/2;
   else
       faceCentersVec(1,:,i) = center - Ztmp*24/2;
       faceCentersVec(2,:,i) = center + Ztmp*24/2;
   end
   
   if sensorPoints
       % Find sensor position by finding the face center nearest another vertex
       d1 = min(sum((V - faceCenters(1,:)).^2, 2));
       d2 = min(sum((V - faceCenters(2,:)).^2, 2));

       if d1 < d2
           Z(i,:) = -Ztmp;
           pos(i,:) = faceCenters(1,:) - 6.5*Z(i,:);
           Y(i,:) = -Ytmp;
       else
           Z(i,:) = Ztmp;
           pos(i,:) = faceCenters(2,:) - 6.5*Z(i,:);
           Y(i,:) = Ytmp;
       end
   else
       Z(i,:) = Ztmp;
       Y(i,:) = Ytmp;
   end
end

%% Manually set position of bottom face if sensors are pointed
if ~sensorPoints
    
    % Set face center 1 as the one nearest the center of the headcast
    % volume
    headcastMiddle = nanmean(nanmean(faceCentersVec, 3), 1);
    
    d1 = squeeze(sqrt(sum((faceCentersVec(1,:,:) - repmat(headcastMiddle, 1, 1, size(sensorListing,1))).^2, 2)));
    d2 = squeeze(sqrt(sum((faceCentersVec(2,:,:) - repmat(headcastMiddle, 1, 1, size(sensorListing,1))).^2, 2)));
    idx = find(d2 < d1);
    
    Z = -Z;
    Z(idx,:) = -Z(idx,:);
    Y = -Y;
    Y(idx,:) = -Y(idx,:);
    pos = squeeze(faceCentersVec(1,:,:))';
    pos(idx,:) = squeeze(faceCentersVec(2,:,idx))';
    pos = pos - 6.5*Z;
end

%% Write to file

A = [1:size(pos,1); pos'; X'; Y'; Z'];
fileID = fopen(fullfile(savePath, 'Headcast_allAxes.txt'),'w');
fprintf(fileID,'%6s %12s %12s %12s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n','SlotNo', 'Position_x', 'Position_y', 'Position_z', 'X_x', 'X_y', 'X_z', 'Y_x', 'Y_y', 'Y_z', 'Z_x', 'Z_y', 'Z_z');
fprintf(fileID,'%6.f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',A);
fclose(fileID);


