function [decimated_headshape] = decimate_headshape(headshape, method)

% Convert to a pointcloud in MATLAB
headshape_pc = pointCloud(headshape.pos);

switch method
    case 'gridaverage'
        try
    if strcmp(headshape.unit,'mm')
         decimated_headshape = pcdownsample(headshape_pc,'gridAverage',20);
    else
        decimated_headshape = pcdownsample(headshape_pc,'gridAverage',2);
    end 
        catch
            disp('Could not determine headshape units');
        end
        
    case 'nonuniform'
        decimated_headshape = pcdownsample(headshape_pc,...
            'nonuniformGridSample',20);
end

decimated_headshape = decimated_headshape.Location;

end











