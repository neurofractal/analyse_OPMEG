function [data_dssp] = ft_dssp_window(cfg,datain)
% Function to apply the DSSP algorithm on overlapping
% windows of electrophysiological data arranged in a Fieldtrip structure
%
% EXAMPLE USEAGE:   data_dssp = ft_dssp_window(cfg,data)
% ...where, cfg is the input structure
% 
%   cfg.sourcemodel      = structure, source model with precomputed leadfields, see FT_PREPARE_LEADFIELD
%   cfg.dssp             = structure with parameters that determine the behavior of the algorithm
%   cfg.dssp.n_space     = 'all', or scalar. Number of dimensions for the
%                          initial spatial projection.
%   cfg.dssp.n_in        = 'all', or scalar. Number of dimensions of the
%                          subspace describing the field inside the ROI.
%   cfg.dssp.n_out       = 'all', or scalar. Number of dimensions of the
%                          subspace describing the field outside the ROI.
%   cfg.dssp.n_intersect = scalar (default = 0.9). Number of dimensions (if
%                          value is an integer>=1), or threshold for the
%                          included eigenvalues (if value<1), determining
%                          the dimensionality of the intersection.
%   cfg.winsize          = size of the window you wish to apply the DSSP
%                          (in seconds, default = 10)
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

% Authors: Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% check the input data
datain = ft_checkdata(datain, 'datatype', {'raw'}); % FIXME how about timelock and freq?

cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'grid',    'sourcemodel'});

% get the options
cfg.trials       = ft_getopt(cfg, 'trials',  'all', 1);
cfg.channel      = ft_getopt(cfg, 'channel', 'all');
cfg.sourcemodel  = ft_getopt(cfg, 'sourcemodel');
cfg.dssp         = ft_getopt(cfg, 'dssp');         % sub-structure to hold the parameters
cfg.dssp.n_space = ft_getopt(cfg.dssp, 'n_space', 'all'); % number of spatial components to retain from the Gram matrix
cfg.dssp.n_in    = ft_getopt(cfg.dssp, 'n_in', 'all');    % dimensionality of the Bin subspace to be used for the computation of the intersection
cfg.dssp.n_out   = ft_getopt(cfg.dssp, 'n_out', 'all');   % dimensionality of the Bout subspace to be used for the computation of the intersection
cfg.dssp.n_intersect = ft_getopt(cfg.dssp, 'n_intersect', 0.9); % dimensionality of the intersection
cfg.output       = ft_getopt(cfg, 'output', 'original');

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});
datain = ft_selectdata(tmpcfg, datain);
% restore the provenance information
%[cfg, datain] = rollback_provenance(cfg, datain);

% match the input data's channels with the labels in the leadfield
sourcemodel = cfg.sourcemodel;
if ~isfield(sourcemodel, 'leadfield')
  ft_error('cfg.sourcemodel needs to contain leadfields');
end
[indx1, indx2] = match_str(datain.label, sourcemodel.label);
if ~isequal(indx1(:),(1:numel(datain.label))')
  ft_error('unsupported mismatch between data channels and leadfields');
end
if islogical(sourcemodel.inside)
  inside = find(sourcemodel.inside);
else
  inside = sourcemodel.inside;
end
for k = inside(:)'
  sourcemodel.leadfield{k} = sourcemodel.leadfield{k}(indx2,:);
end

% compute the Gram-matrix of the supplied forward model
lf = cat(2, sourcemodel.leadfield{:});
G  = lf*lf';

% Calculate window size in terms of number of data points
wsize = cfg.winsize*datain.fsample;
fprintf('Using %d data-points per window\n',wsize);

% Make copy of the data
data_dssp    = datain;

% For every trial
for trial = 1:length(datain.trial)
    % Get data for this trial
    x = datain.trial{trial};
    
    % Create array of zeros for data
    y=zeros(size(x));
    % Create array of zeros for triangular weighting
    a=zeros(1,size(x,2));
    
    % Weighting? I could probably take this out later
    w=ones(size(x));
    if size(w,1)==1; w=repmat(w,1,size(x,1)); end
    
    
    % Start at 0
    offset=0;
    while true
        disp(offset);
        
        start=offset+1;
        stop=min(size(x,2),offset+wsize);
        
        % if not enough valid samples grow window:
        counter=5;
        while any (sum(min(w(start:stop),1))) <wsize
            if counter <= 0 ; break; end
            start=max(1,start-wsize/2);
            stop=min(size(x,2),stop+wsize/2);
            counter=counter-1;
        end
        if rem(stop-start+1,2)==1; stop=stop-1; end
        wsize2=stop-start+1;
        
        % Do DSSP
        [yy, Ae, N, Nspace, Sout, Sin, Sspace, S] = ...
            dssp(x(:,start:stop), G, cfg.dssp.n_in, cfg.dssp.n_out,...
            cfg.dssp.n_space, cfg.dssp.n_intersect);
        
        % triangular weighting
        if start==1
            b=[ones(1,wsize2/2)*wsize2/2, wsize2/2:-1:1];
        elseif stop==size(x,2)
            b=[1:wsize2/2, ones(1,wsize2/2)*wsize2/2];
        else
            b=[1:wsize2/2, wsize2/2:-1:1];
        end
        
        % overlap-add to output
        y(:,start:stop)=y(:,start:stop)+bsxfun(@times,yy,b);
        
        % Add triangular weighting to output
        a(1,start:stop)=a(start:stop)+b;
        
        % Adjust offset parameter
        offset=offset+wsize/2;
        
        % If we have reached the end of the data BREAK 
        if offset>size(x,2)-wsize/2; break; end
    end
    
    % Adjust triangular weighting
    y=bsxfun(@times,y,1./a); 
    % Find any NaN values and convert to 0
    y(isnan(y))=0;
    
    % Put the data back into Fieldtrip structure
    data_dssp.trial{trial} = y;
end
end


