function [headshape_downsampled] = downsample_headshape_FIL(cfg,...
    headshape)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% headshape_downsampled_new: a new function to downsample headshape 
% information for more accurate coregistration. 
%
% Typically you want the headshape information to contain around 200 
% scalp surface points, alongside discernible facial information 
% (eyebrows, eye-sockets and nose)
%
% Designed for data from a Polhemus system or the new Stucture Sensor
%
% Author: Robert Seymour (robert.seymour@mq.edu.au)
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
%
% - path_to_headshape     = path to .hsp file
% - cfg.facial_info       = 'yes' or 'no';
% - cfg.remove_zlim       = number (remove points which lie 
%                           Xmm above nasion on the z-axis (up-down)
% - cfg.method            = method for decimation ('gridaverage' or
%                           'nonuniform'). Gridaverage seems to give more 
%                           consistent results.
%
%%%%%%%%%%%%%%%%
% Fancy Options:
%%%%%%%%%%%%%%%%
%
% cfg.downsample_facial_info        = 'yes' or 'no' : Do you want to 
%                                   downsample the facial information? 
% cfg.downsample_facial_info_amount = amount of facial info downsampling.
%                                   The higher the number the more 
%                                   downsampling (default = 9).
% cfg.facial_info_above_z           = remove points Xmm above the nasion 
%                                   on the z-axis (up-down)
% cfg.facial_info_below_z           = remove points Xmm below the nasion 
%                                   on the z-axis (up-down)
% cfg.facial_info_above_y           = remove points Xmm above the nasion 
%                                   on the y-axis (left-right) 
% cfg.facial_info_above_y           = remove points Xmm below the nasion 
%                                   on the y-axis (left-right) 
% cfg.facial_info_below_x           = remove points Xmm below the nasion 
%                                   on the x-axis (forwards-backwards)
%
%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%
%
% - downsampled_headshape = the downsampled headshape information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default Options
% If not specified include the facial points

facial_info            = ft_getopt(cfg,'facial_info','yes');
if strcmp(facial_info,'no')
    facial_info = [];
end

remove_zlim            = ft_getopt(cfg,'remove_zlim',0);

method                 = ft_getopt(cfg,'method','gridaverage');

facial_info_above_z            = ft_getopt(cfg,'facial_info_above_z',20);
facial_info_below_z            = ft_getopt(cfg,'facial_info_below_z',80);
facial_info_above_y            = ft_getopt(cfg,'facial_info_above_y',70);
facial_info_below_y            = ft_getopt(cfg,'facial_info_above_y',70);
facial_info_below_x            = ft_getopt(cfg,'facial_info_below_x',20);
downsample_facial_info         = ft_getopt(cfg,...
    'downsample_facial_info',[]);

if strcmp(downsample_facial_info,'no')
    downsample_facial_info = [];
end

downsample_facial_info_amount  = ft_getopt(cfg,...
    'downsample_facial_info_amount',9);

%% Headshape
% Convert to mm
headshape = ft_convert_units(headshape,'mm');
% Save a version for later
headshape_orig = headshape;

%% If preserving facial info...
if ~isempty(facial_info)
    
    count_facialpoints = find(headshape.pos(:,3)<facial_info_above_z ...
        & headshape.pos(:,3)>-facial_info_below_z ...
        & headshape.pos(:,1)>facial_info_below_x ...
        & headshape.pos(:,2)<facial_info_above_y ...
        & headshape.pos(:,2)>-facial_info_below_y);
    
    if isempty(count_facialpoints)
        warning('CANNOT FIND ANY FACIAL POINTS - COREG MAY BE INACCURATE');
    else
        facialpoints = headshape.pos(count_facialpoints,:,:);
    end
    % Remove facial points for now
    headshape.pos(count_facialpoints,:) = [];
    
    % If the user specified to downsample the facial info, do this!
    
    if ~isempty(downsample_facial_info) && ~isempty(count_facialpoints)
        
        headshape_pc = pointCloud(facialpoints);
        decimated_facial_info = pcdownsample(headshape_pc,...
            'gridAverage',downsample_facial_info_amount);
        
        fprintf('Downsampling facial info from %d points to %d points\n',...
            length(facialpoints),length(decimated_facial_info.Location));

        facialpoints = decimated_facial_info.Location;
        clear headshape_pc decimated_facial_info
    end
    
    % Plot the facial and head points in separate colours
    figure;
    if isempty(count_facialpoints)
        disp('Not plotting any facial points')
    else
        ft_plot_mesh(facialpoints,'vertexcolor','r','vertexsize',10); hold on;
    end
    ft_plot_headshape(headshape,'vertexcolor','k','vertexsize',10); hold on;
    view([90 0]);
end

%% if the user specified to remove points along the zlim
if ~isempty(remove_zlim)
    
    % Also remove points Xmm above nasion
    points_below_X = find(headshape.pos(:,3) < remove_zlim);
    
    if isempty(points_below_X)
        warning('CANNOT FIND ANY POINTS LESS THAN the zlim specified');
    else
        
        % Remove these points
        headshape.pos(points_below_X,:) = [];
        
        % Plot the mesh again
        figure;
        ft_plot_mesh(headshape.pos,'vertexcolor','r','vertexsize',10); hold on;
        
    end
end

% Create mesh out of headshape downsampled to x points specified in the
% function call
cfg = [];
%cfg.numvertices = 1000;
cfg.method = 'headshape';
cfg.headshape = headshape.pos;
mesh = ft_prepare_mesh(cfg, headshape);

%% Decimate the headshape
try
    [decimated_headshape] = decimate_headshape(headshape, method);
catch
    disp('Incorrect Method Specified');
end

%% Create figure for quality checking
figure; subplot(2,2,1);ft_plot_mesh(mesh,'facecolor','k',...
    'facealpha',0.1,'edgealpha',0); hold on;
ft_plot_mesh(headshape_orig.pos,'vertexcolor','r','vertexsize',2); hold on;
ft_plot_mesh(decimated_headshape,'vertexcolor','b','vertexsize',10); hold on;
view(-180,0);
subplot(2,2,2);ft_plot_mesh(mesh,'facecolor','k',...
    'facealpha',0.1,'edgealpha',0); hold on;
ft_plot_mesh(headshape_orig.pos,'vertexcolor','r','vertexsize',2); hold on;
ft_plot_mesh(decimated_headshape,'vertexcolor','b','vertexsize',10); hold on;
view(0,0);
subplot(2,2,3);ft_plot_mesh(mesh,'facecolor','k',...
    'facealpha',0.1,'edgealpha',0); hold on;
ft_plot_mesh(headshape_orig.pos,'vertexcolor','r','vertexsize',2); hold on;
ft_plot_mesh(decimated_headshape,'vertexcolor','b','vertexsize',10); hold on;
view(90,0);
subplot(2,2,4);ft_plot_mesh(mesh,'facecolor','k',...
    'facealpha',0.1,'edgealpha',0); hold on;
ft_plot_mesh(headshape_orig.pos,'vertexcolor','r','vertexsize',2); hold on;
ft_plot_mesh(decimated_headshape,'vertexcolor','b','vertexsize',10); hold on;
view(-90,0);

print('headshape_quality','-dpng');

% Replace headshape.pos with decimated pos
headshape.pos = decimated_headshape;

% Add the facial points back in (default) or leave out if user specified
% 'no' in function call
if ~isempty(facial_info)
    try
        % Add the facial info back in
        headshape.pos = vertcat(headshape.pos,facialpoints);
    catch
        warning('Cannot add facial info back into headshape');
    end
else
    headshape.pos = headshape.pos;
    warning('Not adding facial points back into headshape');
end

% Plot for quality checking

view_angle = [-180, 0];
figure;

for angle = 1:length(view_angle)
    subplot(1,2,angle)
    ft_plot_headshape(headshape,'vertexcolor','k','vertexsize',12) %plot headshape
    hold on;
    ft_plot_headshape(headshape_orig,'vertexcolor','r','vertexsize',1) %plot headshape
    view(view_angle(angle),10);
end

print('headshape_quality2','-dpng');


% Export filename
headshape_downsampled = headshape;
