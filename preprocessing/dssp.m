%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions for the computation of the projection matrix
% kindly provided by Kensuke, and adjusted a bit by Jan-Mathijs
function [Bclean, Ae, Nee, Nspace, Sout, Sin, Sspace, S] = dssp(B, G, Nin, Nout, Nspace, Nee)

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

[Ae, Nee, Sout, Sin, S] = CSP01(Bin, Bout, Nin, Nout, Nee);
Bclean      = B - (B*Ae)*Ae'; % avoid computation of Ae*Ae'

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
