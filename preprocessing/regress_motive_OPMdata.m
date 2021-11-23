function [dataout] = regress_motive_OPMdata(rawData_MEG,ref,win_size)
  % Function for regressing the Motive data from the (raw) OPM data  
  
    megind          = rawData_MEG.trial{1};
    megres          = zeros(size(megind));
        
    ref_size        = size(ref,1);
    winSize         = size(megind,2);
    
    % add a mean column to the reference regressors
    intercept       = ones(winSize,1);
    
    reference       = ref;
    reference       = [reference ones(size(reference,1),1)];
    
    % Start of Window-Loop
    % Calculate window size in terms of number of data points
    wsize = win_size*rawData_MEG.fsample;
    fprintf('Using %d data-points per window\n',wsize);
    
    % Create array of zeros for data
    meg_ind_synth_grad  = zeros(size(megind));
    % Create array of zeros for triangular weighting
    a                   = zeros(1,size(megind,2));
    
    % Weighting? I could probably take this out later
    w=ones(size(megind));
    if size(w,1)==1; w=repmat(w,1,size(megind,1)); end
    
    % Start at 0
    offset=0;
    ft_progress('init', 'etf', 'Regressing...')
    
    while true
        ft_progress(offset/size(megind,2))
        
        % Calculate start and stop points for this loop
        start=offset+1;
        stop=min(size(megind,2),offset+wsize);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % These values are grown by a factor of 5
        % Is this needed??
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        counter=0;
        while any (sum(min(w(start:stop),1))) <wsize
            if counter <= 0 ; break; end
            start=max(1,start-wsize/2);
            stop=min(size(megind,2),stop+wsize/2);
            counter=counter-1;
        end
        if rem(stop-start+1,2)==1; stop=stop-1; end
        wsize2=stop-start+1;
        
        % Do gradiometry
        beta    = pinv(reference(start:stop,:))*megind(:,start:stop)';
        yy      = (megind(:,start:stop)' - reference(start:stop,:)*beta)';
        
        % triangular weighting (specified via b variable)
        if start==1
            b=[ones(1,wsize2/2)*wsize2/2, wsize2/2:-1:1];
        elseif stop==size(megind,2)
            b=[1:wsize2/2, ones(1,wsize2/2)*wsize2/2];
        else
            b=[1:wsize2/2, wsize2/2:-1:1];
        end
        
        % Add to meg_ind_synth_grad variable outside the loop, weighted by b
        meg_ind_synth_grad(:,start:stop)=meg_ind_synth_grad(:,start:stop)+bsxfun(@times,yy,b);
        
        % Add triangular weighting to variable outside the loop
        a(1,start:stop)=a(start:stop)+b;
        
        % Adjust offset parameter by window size divided by 5
        offset=offset+wsize/5;
        
        % If we have reached the end of the data BREAK
        if offset>size(megind,2)-wsize/5; break; end
    end
    ft_progress('close')
    
    % Adjust for triangular weighting
    meg_ind_synth_grad=bsxfun(@times,meg_ind_synth_grad,1./a);
    % Find any NaN values and convert to 0
    meg_ind_synth_grad(isnan(meg_ind_synth_grad))=0;
    
    dataout            = rawData_MEG;
    dataout.trial{1}   = meg_ind_synth_grad;
end
    