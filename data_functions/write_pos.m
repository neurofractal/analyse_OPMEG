function write_pos(path,data)

% writes out a .pos (polhemus) file with containing digitised scalp and
% fiducial marker data.

% this always seems to be in mm
data = ft_convert_units(data, 'mm');
[p, f, x] = fileparts(path);
filename = fullfile(p, [f '.pos']);
fid = fopen(filename, 'wt');
% write out number of points
fprintf(fid, ' %d\n', length(data.pos));
% write out each point 
for i=1:length(data.pos)
    fprintf(fid, ' %-3d%27.14g%27.14g%27.14g\n', i, data.pos(i,1), data.pos(i,2), data.pos(i,3));
end
if isfield(data,'fid')
    for i=1:numel(data.fid.label)
            fprintf(fid, ' %-s%27.14g%27.14g%27.14g\n', data.fid.label{i}, data.fid.pos(i,1), data.fid.pos(i,2), data.fid.pos(i,3));
    end
end
fclose(fid);