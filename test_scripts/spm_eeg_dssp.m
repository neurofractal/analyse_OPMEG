function Dssp = spm_eeg_dssp(S)

% Robert Seymour, George O'Neill
% $Id

SVNrev = '$Rev: 0000 $';
spm('FnBanner', mfilename,SVNrev);

% Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'D'),                error('Please specify a dataset');  end
if ~isfield(S, 'winsize'),              S.winsize = 10;                 end
if ~isfield(S, 'ssp'),                  S.ssp = 0;                      end
if ~isfield(S, 'dssp'),                 S.dssp = struct;                end
if ~isfield(S.dssp,'n_space'),          S.dssp.n_space = 'all';         end
if ~isfield(S.dssp,'n_in'),             S.dssp.n_in = 'all';            end
if ~isfield(S.dssp,'n_out'),            S.dssp.n_out = 'all';           end
if ~isfield(S.dssp,'n_intersect'),      S.dssp.n_intersect = 0.9;       end
if ~isfield(S,'prefix'),                S.prefix = 'x';                 end


% Load and check if source model specification and forward computation OK
%--------------------------------------------------------------------------

D = spm_eeg_load(S.D);
[~,val] = spm_eeg_inv_check(D);

if ~val
    error('No sourcemodel prepared!')
else
    if ~isfield(D.inv{val},'gainmat')
        error('No gain matrix prepared, please run spm_eeg_lgainmat')
    else
        L = spm_load(fullfile(D.path,D.inv{val}.gainmat));
        % Generate gram matrix of forward solutions
        G = L.G * L.G';
        % Indetify which channels we need (i.e ones with associated fwds)
        chans = match_str(D.chanlabels,L.label);
    end
end

wsize = S.winsize*D.fsample;
fprintf('%-40s: %30s\n','Samples per window',num2str(wsize));


% Load data into memory for faster computation later
%--------------------------------------------------------------------------
Dtmp = D(:,:,:);
Dtmp(chans,:,:) = 0;

ntrials = size(D,3);


% Perform dSSP
%--------------------------------------------------------------------------
spm('FigName','M/EEG dSSP'); spm('Pointer', 'Watch');
spm_progress_bar('Init', ntrials,'Applying dSSP'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
elseif ntrials == 1, Ibar = 1:100;
else Ibar = 1:ntrials; end
fprintf('Perfoming dSSP \n');

for ii = 1:ntrials
    
    % load trial data
    x = D(chans,:,ii);
    % Create array of zeros for data output
    y=zeros(size(x));
    % Create array of zeros for triangular weighting
    a = zeros(1,size(x,2));
    
    % Weighting? I could probably take this out later - RS
    w=ones(size(x));
    if size(w,1)==1; w=repmat(w,1,size(x,1)); end
    
    % Start at 0
    offset=0;
    
    while true
        
        start=offset+1;
        stop=min(size(x,2),offset+wsize);
        
        if rem(stop-start+1,2)==1; stop=stop-1; end
        wsize2=stop-start+1;
        
        % Do DSSP
        yy = ...
            dssp(x(:,start:stop), G, S.dssp.n_in, S.dssp.n_out,...
            S.dssp.n_space, S.dssp.n_intersect, S.ssp);
        
        % triangular weighting
        if start==1
            b=[ones(1,wsize2/2)*wsize2/2, wsize2/2:-1:1];
        elseif stop==size(x,2)
            b=[1:wsize2/2, ones(1,wsize2/2)*wsize2/2];
        else
            b=[1:wsize2/2, wsize2/2:-1:1];
        end
        
        % overlap-add to output
        y(:,start:stop) = y(:,start:stop) + bsxfun(@times,yy,b);
        
        % Add triangular weighting to output
        a(1,start:stop)=a(start:stop)+b;
        
        
        % Adjust offset parameter
        offset=offset+wsize/2;
        
        % If we have reached the end of the data BREAK
        if offset>size(x,2)-wsize/2; break; end
        
        if ntrials == 1
            prog = offset/size(D,2);
            if randi(4,1) == 4
                spm_progress_bar('Set', prog); drawnow;
            end
        end
        
    end
    
    % Adjust triangular weighting
    y = bsxfun(@times,y,1./a);
    % Find any NaN values and convert to 0
    y(isnan(y))=0;
    
    % put data back into array
    Dtmp(chans,:,ii) = y;
    
    if ntrials > 1
        if ismember(ii, Ibar)
            spm_progress_bar('Set', ii); drawnow;
        end
    end
    
end

% Save
%--------------------------------------------------------------------------
Dssp = clone(D,[S.prefix D.fname]);
Dssp(:,:,:) = Dtmp;
save(Dssp);

%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm_progress_bar('Clear');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions for the computation of the projection matrix
% kindly provided by Kensuke, and adjusted a bit by Jan-Mathijs
function [Bclean, Ae, Nee, Nspace, Sout, Sin, Sspace, S] = dssp(B, G, Nin, Nout, Nspace, Nee, ssp_only)

% Nc: number of sensors
% Nt: number of time points
% inputs
% B(Nc,Nt):  interference overlapped sensor data
% G(Nc,Nc): Gram matrix of voxel lead field
% Nout and Nin: dimensions of the two row spaces
% recom_Nspace: recommended value for the dimension of the pseudo-signal subspace
% outputs
% Bclean(Nc,Nt): cleaned sensor data
% Nee: dimension of the intersection
% Nspace: dimension of the pseudo-signal subspace
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%
% The code below is modified by Jan-Mathijs, no functional changes
% merely cosmetics

% eigen decomposition of the Gram matrix, matrix describing the spatial
% components

[U,S]   = eig(G);
Sspace  = abs(diag(S));

[Sspace, iorder] = sort(-Sspace);
Sspace           = -Sspace;
U(:,:)           = U(:,iorder);

if isempty(Nspace)
    ttext = 'enter the spatial dimension: ';
    Nspace    = input(ttext);
elseif ischar(Nspace) && isequal(Nspace, 'interactive')
    figure, plot(log10(Sspace),'-o');
    Nspace = input('enter spatial dimension of the ROI subspace: ');
elseif ischar(Nspace) && isequal(Nspace, 'all')
    Nspace = find(Sspace./Sspace(1)>1e5*eps, 1, 'last');
end
%fprintf('Using %d spatial dimensions\n', Nspace);

% spatial subspace projector
Us   = U(:,1:Nspace);
USUS = Us*Us';

% Bin and Bout creations
Bin  =                  USUS  * B;
Bout = (eye(size(USUS))-USUS) * B;

% create the temporal subspace projector and apply it to the data
%[AeAe, Nee] = CSP01(Bin, Bout, Nout, Nin, Nee);
%Bclean      = B*(eye(size(AeAe))-AeAe);

if ~ssp_only
    Bclean = Bin;
else
    [Ae, Nee, Sout, Sin, S] = CSP01(Bin, Bout, Nin, Nout, Nee);
    Bclean      = B - (B*Ae)*Ae'; % avoid computation of Ae*Ae'
end

function [Ae, Nee, Sout, Sin, S] = CSP01(Bin, Bout, Nin, Nout, Nee)
%
% interference rejection by removing the common temporal subspace of the two subspaces
% K. Sekihara,  March 28, 2012
% Golub and Van Loan, Matrix computations, The Johns Hopkins University Press, 1996
%
%  Nc: number of channels
%  Nt: number of time points
% inputs
%  Bout(1:Nc,1:Nt): interference data
%  Bin(1:Nc,1:Nt): signal plus interference data
%  ypost(1:Nc,1:Nt): denoised data
%  Nout: dimension of the interference subspace
%  Nin: dimension of the signal plus interference subspace
%  Nee: dimension of the intersection of the two subspaces
% outputs
% Ae = matrix from which the projector onto the intersection can
%      be obtained:
% AeAe: projector onto the intersection, which is equal to the
%       interference subspace.
% Nee: dimension of the intersection
%  ------------------------------------------------------------
%  programmed by K. Sekihara,  Signal Analysis Inc.
%  All right reserved by Signal Analysis Inc.
% -------------------------------------------------------------
%

[dum,Sout,Vout] = svd(Bout,'econ');
[dum,Sin, Vin]  = svd(Bin, 'econ');
Sout = diag(Sout);
Sin  = diag(Sin);

if isempty(Nout)
    ttext = 'enter the spatial dimension for the outside field: ';
    Nout  = input(ttext);
elseif ischar(Nout) && isequal(Nout, 'interactive')
    figure, plot(Sout,'-o');
    Nout = input('enter dimension of the outside field: ');
elseif ischar(Nout) && isequal(Nout, 'all')
    Nout = find(Sout./Sout(1)>1e5*eps, 1, 'last');
end
%fprintf('Using %d spatial dimensions for the outside field\n', Nout);

if isempty(Nin)
    ttext = 'enter the spatial dimension for the inside field: ';
    Nin  = input(ttext);
elseif ischar(Nin) && isequal(Nin, 'interactive')
    figure, plot(log10(Sin),'-o');
    Nin = input('enter dimension of the inside field: ');
elseif ischar(Nin) && isequal(Nin, 'all')
    Nin = find(Sin./Sin(1)>1e5*eps, 1, 'last');
end
%fprintf('Using %d spatial dimensions for the inside field\n', Nin);

Qout = Vout(:,1:Nout);
Qin  = Vin(:, 1:Nin);

C     = Qin'*Qout;
[U,S] = svd(C);
S     = diag(S);
if (ischar(Nee) && strcmp(Nee, 'auto'))
    ft_error('automatic determination of intersection dimension is not supported');
elseif ischar(Nee) && strcmp(Nee, 'interactive')
    figure, plot(S,'-o');
    Nee  = input('enter dimension of the intersection: ');
elseif Nee<1
    % treat a numeric value < 1 as a threshold
    Nee = find(S>Nee,1,'last');
    if isempty(Nee), Nee = 1; end
end
%fprintf('Using %d dimensions for the interaction\n', Nee);

Ae   = Qin*U;
Ae   = Ae(:,1:Nee);
%AeAe = Ae*Ae';