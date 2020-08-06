function spect = spectrogram_emd(F,A,dim)
% This function produces three dimensional spectrum matrix to plot the
% instantaneous amplitude versus time and instantaneous frequency, that is,
% amplitude contours on the time-frequency plane.
%
% INPUT
% F : [N X M] instantaneous frequency matrix
%     The value in the matrix varies between 0 and 0.5.
%     (N : the number of IMF, M : data length)
% A : [N X M] instantaneous amplitude matrix
%     (N : the number of IMF, M : data length)
% dim : dimension of frequency domain
%       (This defines the width of frequency bin. Normally value of dim=1000 is used.)
%
% OUTPUT
% spect : [dim X M] Hilbert-Huang spectrum (HHS)
%         (dim : the dimension of frequency domain, defined by the input `dim',
%          M : data length)


%F_idx= (Fs*.5)/dim:(Fs*.5)/dim:Fs*.5;
%tmp_row=[];tmp_col=[];
[row,col] = size(F);
%t_idx = 1/Fs:1/Fs:col/Fs;
spect = zeros(dim,col);

g = 0:0.5/(dim):.5-0.5/(dim);
h = repmat(g,[row,1]);

for j=1:col
    F_rep = repmat(F(:,j),[1 dim]); % repeat F, 'dim'-times (columns)
    
    % Perform find simultaneously of (F_rep and h)
    [I,J] = find(F_rep>=h & F_rep<(h+1/(2*dim))); % I = row, J = column
    
    % Update spect
    leng = length(I);
    for k=1:leng
        spect(J(k),j) = spect(J(k),j)+A(I(k),j);
    end
end
