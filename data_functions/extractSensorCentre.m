function sensorInfo = extractSensorCentre(cfg)
% Extracts the orientations (rad and tan, sensor space)and the physical
% centre of the sensor. This is an alternative to the centre of the cell.
% 
% This code is quite specific to the STL files for the Gen 2 OPMs, but
% could be adapted for others.
% 
% Inputs. All are required at this stage, but can input empty. 
% cfg.folder        = 'string', folder containing only the STL for sensors.
% cfg.slotori       = 'rad' or 'tan', the sensor-space radial to the head.
% cfg.plot          = 'yes' or 'no'
% cfg.correct       = 'yes' or 'no' : do you want to manually flip the tan ori?
% cfg.scalp         = 'string', path to STL file of scalp in the same coordinate system.
% cfg.outputfolder  = 'string', filename with directory for output.
% cfg.sensortype    = 'old' or 'new', where old is the plain rectangle with
%                       cone and new is the accurate model of the sensor.
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
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Edited from Stephanie Mellor's (09/02/2020) code by Nicholas Alexander 
% (12/03/2020) and Robert Seymour (16/03/2020)
%__________________________________________________________________________

%% Example Usage
% cfg               = [];
% cfg.folder        = '/Volumes/Robert T5/OPM_data/gareth_scanner_cast/indiv_stl';
% cfg.slotori       = 'tan';
% cfg.plot          = 'yes';
% cfg.correct       = 'yes';
% cfg.outputfolder  = '/Volumes/Robert T5/OPM_data/gareth_scanner_cast/hingecast';
% cfg.scalp         = '/Volumes/Robert T5/OPM_data/gareth_scanner_cast/Head.stl';
% sensorInfo        = extractSensorCentre(cfg)


%% Start of function
% Input file path and find all STL files. 

close all force

STLdir              = strcat(cfg.folder);

sensorListing       = dir([STLdir '\*.stl']);

if isempty(sensorListing)
    sensorListing       = dir([STLdir '/*.stl']);
end


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
filename            = cell(size(sensorListing, 1), 1);

for sensorIdx = 1:size(sensorListing, 1)
    
    fprintf('Extracting information from .stl file %3d  of %3d\n',...
        sensorIdx,size(sensorListing, 1));
    
    % Specifiy and read in sensor STLs one at a time. 
    sensorFilename              = fullfile(sensorListing(sensorIdx).folder,...
        sensorListing(sensorIdx).name);
    [faces, verts]              = stlread(sensorFilename);
    
    filename{sensorIdx}         = sensorFilename;
    
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
    
    
    switch cfg.sensortype
        case 'old'
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

            % Remove the hypotenuses.
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
        case 'new'
            % There is a cone on the top and bottom of the sensor. 
            % The cone vertices can be found by looking for the top two
            % most referenced vertices within the faces of the model. 
            facesList                   = reshape(redFaces,[numel(redFaces),1]);
            [n,bin]                     = hist(facesList,unique(facesList)); %#ok<HIST>
            [~,idx]                     = sort(-n);
            count                       = n(idx); % count instances
            value                       = bin(idx); % corresponding values
            
            % Fortunately, the most common one is always the bottom and 
            % second most is always the top (98 and 99)
            bottomPoint                 = [redVerts(value(1),1),redVerts(value(1),2),redVerts(value(1),3)];
            topPoint                    = [redVerts(value(2),1),redVerts(value(2),2),redVerts(value(2),3)];
            centrePoint(sensorIdx,1:3)  = mean([bottomPoint; topPoint]);
%             
%             % Debug plotting
%             hold on
%             patch('Faces',redFaces,'Vertices',redVerts,'FaceColor','red');
%             scatter3(bottomPoint(1),bottomPoint(2),bottomPoint(3));
%             scatter3(topPoint(1),topPoint(2),topPoint(3));
%             scatter3(centrePoint(1),centrePoint(2),centrePoint(3));
%             hold off
            
            % The shape has internal vertices but we need the outer 
            % boundary.
            boundFaces              = convhull(redVerts,'simplify',true);            
%             
%             % Debug plots.
%             hold on
%             patch('Faces',boundFaces,'Vertices',redVerts,'FaceColor','red');
%             scatter3(redVerts(:,1),redVerts(:,2),redVerts(:,3));
%             
            % Remove unused vertices.
            uniqueBoundFaces        = unique(boundFaces);
            boundVerts              = [];
            newBoundFaces           = zeros(size(boundFaces));
            boundVerts              = zeros(length(uniqueBoundFaces),3);
            
            for vertIdxIdx = 1:length(uniqueBoundFaces)
                boundVerts(vertIdxIdx,1:3) = redVerts(uniqueBoundFaces(vertIdxIdx),1:3);
                
                % Replace all references to the redVerts idx with the
                % boundVerts idx. 
                [xIdx, yIdx] = find(boundFaces == uniqueBoundFaces(vertIdxIdx));
                for ii = 1:length(xIdx)
                    newBoundFaces(xIdx(ii),yIdx(ii)) = vertIdxIdx;
                end
            end
            
%             % debug plot
%             patch('Faces',newBoundFaces,'Vertices',boundVerts,'FaceColor','red');

            % Fit cuboid to the boundary points.
            [simpleVerts, cuboidParameters, ~, ~] = CuboidRANSAC(boundVerts);   
            
            simpleFaces     = convhull(simpleVerts,'simplify',true);
            
            % Take the simple cuboid dimensions for later
            edgeDists       = round(cuboidParameters(4:6),4)';
            
%             % debug plot
%             hold on
%             patch('Faces',simpleFaces,'Vertices',simpleVerts,'FaceColor','red');
%             scatter3(redVerts(:,1),redVerts(:,2),redVerts(:,3));
            


            % Find groups of the long sides.
            longestSide             = max(edgeDists);
            edgeDists(edgeDists == longestSide) = [];


            % Take a vertex and find the connected vertices. 
            firstGroupCornerVerts   = zeros(length(simpleVerts(:,1))/2,length(simpleVerts(1,:)));
            remCornerVerts          = simpleVerts;
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
                    dist    = round(dist,4);
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
            
%             % Debug plot
%             hold on
%             scatter3(firstGroupCornerVerts(:,1),firstGroupCornerVerts(:,2),firstGroupCornerVerts(:,3));
%             scatter3(secondGroupCornerVerts(:,1),secondGroupCornerVerts(:,2),secondGroupCornerVerts(:,3));
%             patch('Faces',redFaces,'Vertices',redVerts,'FaceColor','red');
            
            % Take the mean position of each grouping.
            firstGroupCentreVert    = mean(firstGroupCornerVerts,1);
            secondGroupCentreVert   = mean(secondGroupCornerVerts,1);

            % Which is closer to the centre of the sensor?
            firstGroupCentreDist    = pdist([firstGroupCentreVert;centrePoint(sensorIdx,1:3)]);
            secondGroupCentreDist   = pdist([secondGroupCentreVert;centrePoint(sensorIdx,1:3)]);

            if (firstGroupCentreDist < secondGroupCentreDist)
                flatEndVerts        = firstGroupCornerVerts;
                flatEndCentreVert   = firstGroupCentreVert;
                cableEndVerts       = secondGroupCornerVerts;
                cableEndCentreVert  = secondGroupCentreVert;
            elseif (firstGroupCentreDist > secondGroupCentreDist)
                flatEndVerts        = secondGroupCornerVerts;
                flatEndCentreVert   = secondGroupCentreVert;
                cableEndVerts       = firstGroupCornerVerts;
                cableEndCentreVert  = firstGroupCentreVert;
            end   
            
            % Very confusing, but if the sensors are in 'tan' orientation slots,
            % then that means the tan of the sensor is aligned with the rad of the
            % head. Need to avvound for that. 
            if strcmp(cfg.slotori,'rad')
                radOri(sensorIdx,1:3)       = cableEndCentreVert - flatEndCentreVert;
                tanOri(sensorIdx,1:3)       = topPoint - bottomPoint;
            elseif strcmp(cfg.slotori,'tan')
                radOri(sensorIdx,1:3)       = topPoint - bottomPoint;
                tanOri(sensorIdx,1:3)       = cableEndCentreVert - flatEndCentreVert;
            end
            
            % Find the physical centre instead
            centrePoint(sensorIdx,1:3) = mean([firstGroupCornerVerts(:,:); secondGroupCornerVerts(:,:)]);
            
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
if strcmp(cfg.plot,'yes')
    figure(1);
    hold on;
    grid off;
    set(gca,'visible','off')
    quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3),...
        radOri(:,1), radOri(:,2), radOri(:,3));
    quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3),...
        tanOri(:,1), tanOri(:,2), tanOri(:,3));
    scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
    text(centrePoint(:,1)-1.5*radOri(:,1),centrePoint(:,2)-1.5*radOri(:,2),...
        centrePoint(:,3)-1.5*radOri(:,3), cellstr(num2str((1:size(centrePoint,1))')), 'color', 'g');
    daspect([1 1 1])
    for sensorIdx = 1:size(sensorListing, 1)
        patch('Faces',sensorFaces{sensorIdx},'Vertices',sensorVerts{sensorIdx},'FaceAlpha',.3,'EdgeAlpha',0);
    end
    
    if ~isempty(cfg.scalp)
        [scalpFaces, scalpVerts]    = stlread(cfg.scalp);
        cindex                      = scalpVerts(:,3);
        patch('Faces',scalpFaces,'Vertices',scalpVerts,'FaceVertexCData',cindex,'FaceColor','interp','EdgeAlpha',0);
        colormap(copper)
    end
    hold off
    % Ask if they are done with it.
    input('Press any key to continue (closes figure)\n')
    close all % Change this later. 
end

% For determining whether to flip the orientation or not, plot one sensor
% onto the scalp at a time.
if strcmp(cfg.correct,'yes') && strcmp(cfg.sensortype,'old')
    disp('Plotting sensor pos and ori one at a time. Please manually correct...');
    repeat = 1;
    while repeat
        % If user specified a scalp,load this now
        if ~isempty(cfg.scalp)
                [scalpFaces, scalpVerts]  = stlread(cfg.scalp);
        end
        
        % Create Figure
        S.f = figure; 
        
        % For every sensor
        for sensorIdx = 1:size(sensorListing, 1)
            % Plot Scalp mesh if specified
            if ~isempty(cfg.scalp)
                S.h = patch('Faces',scalpFaces,'Vertices',scalpVerts,'FaceVertexCData',cindex,'FaceColor','interp','EdgeAlpha',0);
                colormap(copper)
            end
            hold on;
            %create two pushbttons
            S.pb = uicontrol('style','push',...
                'units','pix',...
                'position',[370 10 180 40],...
                'fontsize',14,...
                'Tag','flip_button',...
                'string','Flip = YES',...
                'UserData',struct('flip',0),...
                'callback',@pb_call);
            
            S.pb = uicontrol('style','push',...
                'units','pix',...
                'position',[10 10 180 40],...
                'fontsize',14,...
                'UserData',struct('flip',0),...
                'string','Flip = NO',...
                'callback',@pb_call2);
            
            S.h = patch('Faces',sensorFaces{sensorIdx},'Vertices',...
                sensorVerts{sensorIdx},'FaceAlpha',1,'EdgeAlpha',...
                0,'FaceColor',[0 0 0]);
            
            hold on;
            
            S.h = quiver3(centrePoint(sensorIdx,1), centrePoint(sensorIdx,2),...
                centrePoint(sensorIdx,3), tanOri(sensorIdx,1).*50,...
                tanOri(sensorIdx,2).*50, tanOri(sensorIdx,3).*50,'color',...
                [0 0 1],'LineWidth',3); 
            
            hold on;
            
            S.h = scatter3(centrePoint(sensorIdx,1), centrePoint(sensorIdx,2),...
                centrePoint(sensorIdx,3));
           
            % Remove axes
            set(gca,'visible','off')
            
            % Set view
            view(-radOri(sensorIdx,:));
            
            % Wait for user input
            uiwait(S.f)
            h = findobj('Tag','flip_button');
            
            % If user pressed 'FLIP' then flip the tan orientation
            if h.UserData.flip == 1
                tanOri(sensorIdx,1:3)   = -tanOri(sensorIdx,1:3);
                disp([num2str(sensorIdx) ' FLIPPED']);
            else
                disp([num2str(sensorIdx) ' NOT FLIPPED']);
            end
            
            % Clear the Figure for the next sensor
            clf(S.f);

        end
        
        % Show a plot of the current orientations.
        figure(2);
        hold on;
        grid off;
        set(gca,'visible','off')
        quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3),...
            radOri(:,1), radOri(:,2), radOri(:,3));
        quiver3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3),...
            tanOri(:,1), tanOri(:,2), tanOri(:,3));
        scatter3(centrePoint(:,1), centrePoint(:,2), centrePoint(:,3))
        text(centrePoint(:,1)-1.5*radOri(:,1),centrePoint(:,2)-1.5*radOri(:,2),...
            centrePoint(:,3)-1.5*radOri(:,3), cellstr(num2str((1:size(centrePoint,1))')), 'color', 'g');
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
        repeat    = input('Would you like to go through them again? 1 for yes, 0 for no');
        close all
    end
end

% Put the main info into one output.
filename    = [filename;filename];
slot        = [(1:length(sensorListing))';(1:length(sensorListing))'];
Px          = [centrePoint(:,1);centrePoint(:,1)];
Py          = [centrePoint(:,2);centrePoint(:,2)];
Pz          = [centrePoint(:,3);centrePoint(:,3)];

% Create list describing whether orientation corresponds to 
% TAN or RAD OPM sensors. Where sensors are placed tan (e.g. for hingecast)
% the sensitive axis will need to be flipped
corresponding_sens = cell(length(sensorListing)*2,1);
if strcmp(cfg.slotori,'tan')
    corresponding_sens(1:length(sensorListing))      = {'RAD'};
    corresponding_sens(length(sensorListing)+1:end)  = {'TAN'};
elseif strcmp(cfg.slotori,'rad')
    corresponding_sens(1:length(sensorListing))      = {'TAN'};
    corresponding_sens(length(sensorListing)+1:end)  = {'RAD'};
else
    corresponding_sens(1:end)                       = {'UNKNOWN'};
end


Ox          = [radOri(:,1);tanOri(:,1)];
Oy          = [radOri(:,2);tanOri(:,2)];
Oz          = [radOri(:,3);tanOri(:,3)];

outputTable = table(filename,slot,corresponding_sens,Px,Py,Pz,Ox,Oy,Oz);

if ~isempty(cfg.outputfolder)
    try
    disp('Writing .tsv file...');
    writetable(outputTable,strcat(cfg.outputfolder,'\positions_centre.tsv'),'Delimiter',...
        'tab','QuoteStrings',false,'FileType', 'text');
    catch
       warning(['Was not able to write the .tsv file... ',...
           'Did you create .tsv file with the same name before?']);
    end
end

sensorInfo      = outputTable;

%% Subfunctions for buttons
function pb_call2(varargin)
%     disp('NOT Flipped');
    uiresume;
    return
end

function pb_call(hObject,~)
    hObject.UserData = struct('flip',1);
%     disp('Flipped');
    uiresume;
    return
end
end

