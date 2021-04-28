function [source] = get_source_pow(data_clean,sourceall,toi)
%
% get_source_pow is a function to extract power values from the output of
% ft_sourceanalysis between specific times of interest (toi). The function
% will average power.
%
% Author: Robert Seymour (May 2019) robert.seymour@mq.edu.au
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
%
% - data_clean              = the clean data used for your source analysis
% - sourceall               = output of ft_sourceanalysis for the entire
%                           time-range of interest (ie. baseline and 
%                           trial period)
% - toi                     = times of interest
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs (saved to dir_name):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - source                  = Fieldtrip source structure with modified  
%                           avg.power field reflecting average "power" 
%                           across the toi

% EXAMPLE FUNCTION CALL:
% [sourceN1] = get_source_pow(data_clean,sourceall,[0.1 0.2])



fprintf('Getting source level power from %.3fs to %.3fs\n',toi(1),toi(2)) 

% average across time the dipole moments within the toi
ind    = find(sourceall.time >=toi(1) & sourceall.time <=toi(2));
tmpmom = sourceall.avg.mom(sourceall.inside);
mom    = sourceall.avg.pow(sourceall.inside);
for ii = 1:length(tmpmom)
    mom(ii) = mean(abs(tmpmom{ii}(ind)));
end

% Create new source structure with 'pow' containing the mean amplitude across the
% time-window used for the channel-level covariance computation
source = sourceall;
source.avg.pow = zeros(length(source.pos),1);
source.avg.pow(source.inside) = abs(mom);

try
source.cfg = rmfield(source.cfg,...
    {'headmodel' 'callinfo'}); % this is removed because it takes up a lot of memory
source = rmfield(source,...
    {'time'}); % this is removed because it takes up a lot of memory

catch
    disp('cfg field has already been removed');
end

end