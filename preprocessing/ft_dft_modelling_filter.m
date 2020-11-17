function [filt, filtFreq] = ft_dft_modelling_filter(cfg, data)
% Function to find peaks in spectrum, model their shape and remove them 
% before transforming back into time domain. 
% 
% This function takes inspiration from the spectrum interpolation method 
% (Leske & Dalal, 2019, NeuroImage 189, 10.1016/j.neuroimage.2019.01.026).
% The aim of this function is similar, but intends to automatically detect
% peaks and their distributions so as to minimise the data removed. 
% 
% Overview of function
% 1)   Data transformed into the frequency domain via using Welch method.
% 2)   Peaks in the specified frequency range are identified and some
%       properties extracted about them.
% 3)   Either independently or as a pair, estimates are made about the
%       distribution of peaks. 
% 4)   For each channel, a better fit is made with some channel-wide
%       constraints (centre and max width (g) of peak)
% 5)   Depending on user input, these peak fits are either (1) subtracted 
%       from the original fft power output or (2) used to determine the
%       edges of each peak and to replace data within those edges with a
%       fitted slope. 
% 4)   The signal is transformed back into the time domain via inverse DFT
%       (iDFT). 
%
% Use as
%   [filt, filtFreq] = ft_preproc_dft_remove_gause(cfg, data)
%   where
%   data                    = Fieldtrip data structure with continuous
%                               data. Usually raw data.
%   cfg                     = Structure containing the following:
%   cfg.foi                 = [start end], freq range to remove peaks 
%                               between e.g. [40 60];
%   cfg.Neighwidth          = width (in Hz) of peaks to evaluate, e.g. 2;
%   cfg.minPeakDistance     = minimum distance in Hz between peaks.
%   cfg.method              = 'leaveSlope' or 'removePeak'. Whether to
%                               subtract the modelled noise, potentially 
%                               leaving signal behind, or reduce power 
%                               uniformly. 
% 	cfg.peakShape           = 'Gaussian' or 'Lorentzian'. Determines the
%                               distribution used to model the peak. 
%                               Lorentzian usually works best.
%   cfg.independentPeaks    = 'yes' or 'no'. Whether to model peaks 
%                               independently or together. Note, at the
%                               moment only two peaks may be modelled
%                               together. 
%   cfg.log                 = 'yes' or 'no'. Whether to log the PSD before
%                               fitting. May help with visualisation. 
%   cfg.strength            = multiplier of g, the variable that determines 
%                             distribution width. Higher strength will
%                             remove more of the peak, but increase the
%                             width of data lost.
% 
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% 
% Authors: Nicholas Alexander      (n.alexander@ucl.ac.uk) 
%          Oliver Alexander
%__________________________________________________________________________


%% Extract info from input structure
timeSeries          = data.trial{1};
nsamples            = length(timeSeries);
samplingFreq        = data.fsample;
fftFreq             = samplingFreq * linspace(0, 1, nsamples);
fftFreq             = fftFreq(1:nsamples/2+1);

% Run the fft
fftData             = fft(timeSeries,nsamples,2);

% To get a better estimation, use the Welch method
welchPow            = zeros(length(data.label),length(fftFreq));
for chanIdx = 1:length(data.label)
    if chanIdx == 1
        [welchTmp, welchFreq]       = pwelch(timeSeries(chanIdx,:),[],0,[],samplingFreq);
        welchPow(chanIdx,:)         = interp1(welchFreq,welchTmp,fftFreq);
    else
        welchTmp                    = pwelch(timeSeries(chanIdx,:),[],0,[],samplingFreq);
        welchPow(chanIdx,:)         = interp1(welchFreq,welchTmp,fftFreq);
    end
end

% Put Welch in the same scale as the fft
window          = hanning(nsamples);
scaledWelchPow  = sqrt(welchPow*samplingFreq*(window'*window));
avgWelchPow     = median(scaledWelchPow,1);

% Tidy up
clear timeSeries samplingFreq nsamples 

% Select pow data for the foi.
foiIdx            	= ~(fftFreq < cfg.foi(1) | fftFreq > cfg.foi(2));
if strcmp(cfg.log,'yes')
    avgFftPowFOI           = log(avgWelchPow(foiIdx));
elseif strcmp(cfg.log,'no')
    avgFftPowFOI           = avgWelchPow(foiIdx);
end
freqFOI             = fftFreq(foiIdx);

% The user will be asked if they are satisfied with peak identification.
firstRun            = true;
peakCount           = 0;
while firstRun || changePeaks
    % Get a staring point for peak prominence in the first run.
    if firstRun
        if strcmp(cfg.independentPeaks,'yes')
            [~,~,~,allPeakProm]   = findpeaks(avgFftPowFOI,freqFOI,'MinPeakDistance',cfg.minPeakDistance);
        elseif strcmp(cfg.independentPeaks,'no')
            [~,~,~,allPeakProm]   = findpeaks(avgFftPowFOI,freqFOI);
        end

        if isempty(allPeakProm)
            error('No peaks detected');
        else
            allPeakProm         = sort(allPeakProm,'descend');
            peakCount           = 1;
            minPeakProminence   = allPeakProm(peakCount);
        end
        
        % Prepare figure on first run
        figure;
    end
    
    % Identify peaks based on prominence. Extract halfheight width.
    if strcmp(cfg.independentPeaks,'yes')
        [~,peakFreq,peakWidth,peakProminence]  = findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'MinPeakDistance',cfg.minPeakDistance,'WidthReference','halfprom');
    elseif strcmp(cfg.independentPeaks,'no')
        [~,peakFreq,peakWidth,peakProminence]  = findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'WidthReference','halfprom');
    end

    
    % Show user the peaks
    if strcmp(cfg.independentPeaks,'yes')
        findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'Annotate','extents','MinPeakDistance',cfg.minPeakDistance,'WidthReference','halfprom');
    elseif strcmp(cfg.independentPeaks,'no')
        findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'Annotate','extents','WidthReference','halfprom');
    end
    
    % If the peaks are independent then the centre of the width should be a
    % better estimate of the peak centre.
    if strcmp(cfg.independentPeaks,'yes')
       % Workaround because findpeaks does not output this info.
       currentAxes          = gca;
       currentLines         = currentAxes.Children;
       widthLines           = currentLines(1).XData';
       widthLines           = widthLines(~isnan(widthLines));
       widthLines           = transpose(reshape(widthLines,2,[]));
       peakFreq             = mean(widthLines,2);
    end
    
    % Ask if the user wants to change the number of peaks
    diffPeaks           = 3;
    while ~any(eq(diffPeaks,[0 1 2]))
        diffPeaks = input('To continue, enter 0\nFor less peaks, enter 1\nFor more peaks, enter 2:\n');
    end

    % Set minPeakProminence based on user input
    if diffPeaks == 1
        if minPeakProminence == max(allPeakProm)
            disp('Already at minimum peaks. Enter either 2 or 0');
            changePeaks     = true;
        else
            peakCount           = peakCount - 1;
            minPeakProminence   = allPeakProm(peakCount);
            changePeaks     = true;
        end
    elseif diffPeaks == 2
        if peakCount == length(allPeakProm)
            disp('No more peaks detected. Consider changing FOI')
            changePeaks     = true;
        else
            peakCount           = peakCount + 1;
            minPeakProminence   = allPeakProm(peakCount);
            changePeaks     = true;
        end
    else
        changePeaks = false;
        close(gcf);
        
        % Sort peaks by prominence
        [peakProminence,peakSortIdx]    = sort(peakProminence,'descend');
        peakFreq                        = peakFreq(peakSortIdx);
        peakWidth                       = peakWidth(peakSortIdx);
    end

    if firstRun
        firstRun = false;
    end
end


% Tidy up
clear firstRun changePeaks minPeakProminence diffPeaks noPeaks...
    initialPeakPromList freqFOI foiIdx smoothPowFOI...
    allPeakProm peakSortIdx peakCount avgFftPowFOI freqFOI


%% Define functions for fitting
% Define function for Gaussian/Lorentzian with and without
% slopes.
justSlope       = @(b, c, x)...
                    (b.*x) + c;
switch cfg.peakShape
    case 'Gaussian'
        peakWithSlope   = @(A, x0, g, b, c, x)...
                            (A*exp(-2*((x-x0)/g).^2)) + (b*x) + c;
        justPeak        = @(A, x0, g, x)...
                            (A*exp(-2*((x-x0)/g).^2));
       
    case 'Lorentzian'
        peakWithSlope   = @(A, x0, g, b, c, x)...
                            ((A .* (g.^2))./((g.^2) + ((x - x0).^2))) + (b.*x) + c;
        justPeak        = @(A, x0, g, x)...
                            ((A .* (g.^2))./((g.^2) + ((x - x0).^2)));
                        
        twoPeaksWithSlope = @(A1, x01, g1, A2, x02, g2, b, c, x)...
                                    ((A1 .* (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2))) + (b.*x) + c;

        justTwoPeaks      = @(A1, x01, g1, A2, x02, g2, x)...
                                    ((A1 * (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2)));

end

%% Remove Gaussian component from amplitude in spectrum
switch cfg.independentPeaks
    case 'yes'
        filtFreq        = zeros(length(peakFreq),2);
        for peakIdx = 1:length(peakFreq)
            %% Find the neighbourhood data around the peak
            neighbourLowerBound         = peakFreq(peakIdx) - (cfg.Neighwidth);
            neighbourUpperBound         = peakFreq(peakIdx) + (cfg.Neighwidth);
            neighbourFreqIndicesBound   = nearest(fftFreq,[neighbourLowerBound,neighbourUpperBound]);
            neighbourFreqIndices        = neighbourFreqIndicesBound(1):neighbourFreqIndicesBound(end);
            
            if strcmp(cfg.log,'yes')
                % Get the chan avg data for those bounds.
                neighbourData               = log(scaledWelchPow(:,neighbourFreqIndices));
%                 fftNeighbourData            = log(fftPow(:,neighbourFreqIndices));
            elseif strcmp(cfg.log,'no')
                % Get the chan avg data for those bounds.
                neighbourData               = scaledWelchPow(:,neighbourFreqIndices);
%                 fftNeighbourData            = fftPow(:,neighbourFreqIndices);
            end

            %% First fit the peak with slope to each channel independently
            % Guesses for fit (peak level)
            % Amplitude
            A               = peakProminence(peakIdx);
            
            % Centre
            x0              = nearest(fftFreq,peakFreq(peakIdx));
            
            % Gamma (assumes fftFreq starts at 0)
            g               = nearest(fftFreq,peakWidth(peakIdx));
            
            
            % Get a better estimate of g and x0
            for chanIdx = 1:length(fftData(:,1))
                % Guesses for fit (channel level)
                neighbourChanData   = neighbourData(chanIdx,:);
                % Slope
                quarterLength   = round(length(neighbourChanData)/4);
                endVals         = mean(neighbourChanData((end - quarterLength):end));
                startVals       = mean(neighbourChanData(1:1+quarterLength));
                eighthLength    = round(quarterLength/2);
                b               = (endVals - startVals) / (neighbourFreqIndices(end - eighthLength) - neighbourFreqIndices(eighthLength));
                
                % Constant - What is c equal to for the slope go through the middle?
                middleValue     = mean([endVals,startVals]);
                c               = middleValue - (b * x0);

                % Fit a peak with slope to the channel data.
                fittedModel     = fit(neighbourFreqIndices', neighbourChanData', peakWithSlope,...
                                    'StartPoint', [A, x0, g, b, c],...
                                    'Lower',[-inf, x0-g, -inf, -inf, -inf],...
                                    'Upper',[inf, x0+g, inf, inf, inf]);

                fitCoeffValue   = coeffvalues(fittedModel);
                fitCoeffNames   = coeffnames(fittedModel);
                
                % Dynamic naming seems to be the best way...
                for coeffIdx = 1:length(fitCoeffNames)
                    fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
                end
            end
            
            % Set the x0 for remaining fits
            x0      = median(fittedStruct.x0);
            
            
            % It might be worth using the maxG as a limit for the fitted
            % model. 
            medianG        = median(abs(fittedStruct.g));
            iqrG            = iqr(abs(fittedStruct.g));
            maxG            = medianG + 3*iqrG;
            minG            = abs(medianG - 3*iqrG);
            % Fit again with a contraint on g
            for chanIdx = 1:length(fftData(:,1))
                % Guesses for fit (channel level)
                neighbourChanData   = neighbourData(chanIdx,:);
                % Slope
                quarterLength   = round(length(neighbourChanData)/4);
                endVals         = mean(neighbourChanData((end - quarterLength):end));
                startVals       = mean(neighbourChanData(1:1+quarterLength));
                eighthLength    = round(quarterLength/2);
                b               = (endVals - startVals) / (neighbourFreqIndices(end - eighthLength) - neighbourFreqIndices(eighthLength));
                
                % Constant - What is c equal to for the slope go through the middle?
                middleValue     = mean([endVals,startVals]);
                c               = middleValue - (b * x0);

                % Fit a peak with slope to the channel data.
                fittedModel     = fit(neighbourFreqIndices', neighbourChanData', peakWithSlope,...
                                    'StartPoint', [A, x0, medianG, b, c],...
                                    'Lower',[-inf, x0, minG, -inf, -inf],...
                                    'Upper',[inf, x0, maxG, inf, inf]);

                fitCoeffValue   = coeffvalues(fittedModel);
                fitCoeffNames   = coeffnames(fittedModel);
                
                % Dynamic naming seems to be the best way...
                for coeffIdx = 1:length(fitCoeffNames)
                    fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
                end
            end
            
            maxG            = max(abs(fittedStruct.g));
            
            % Get the width of cfg.strength gamma of fitted lorentzian.
            tmpStartIdx             = nearest(neighbourFreqIndices,(x0 - (cfg.strength * maxG)));
            tmpEndIdx               = nearest(neighbourFreqIndices,(x0 + (cfg.strength * maxG)));

            indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
            
            % Save the freq limits that have been modified. 
            filtFreq(peakIdx,1:2)   = [fftFreq(indicesToReplace(1)),fftFreq(indicesToReplace(end))];

            if length(indicesToReplace) > length(neighbourData(chanIdx,:))
                ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
            end
            
            for chanIdx = 1:length(fftData(:,1))
                switch cfg.method
                    case 'removePeak'
                        % Get just the peak
                        peakOnly                = justPeak(fittedStruct.A(chanIdx),fittedStruct.x0(chanIdx),fittedStruct.g(chanIdx),indicesToReplace);

                        % Remove it from the original data.
                        replacementData         = neighbourData(chanIdx,tmpStartIdx:tmpEndIdx) - peakOnly;
                        
                    case 'removePeak2'
                        % Get the peak on a slope
                        fittedModelData         = peakWithSlope(fittedStruct.A(chanIdx),fittedStruct.x0(chanIdx),fittedStruct.g(chanIdx),fittedStruct.b(chanIdx),fittedStruct.c(chanIdx),indicesToReplace);
                        
                        % Remove it from the original data.
                        replacementData         = neighbourData(chanIdx,tmpStartIdx:tmpEndIdx) - fittedModelData;

                    case 'leaveSlope'
                        % Get the peak on a slope
                        fittedModelData         = peakWithSlope(fittedStruct.A(chanIdx),fittedStruct.x0(chanIdx),fittedStruct.g(chanIdx),fittedStruct.b(chanIdx),fittedStruct.c(chanIdx),neighbourFreqIndices);
                        tmpStartVal             = fittedModelData(tmpStartIdx);
                        tmpEndVal               = fittedModelData(tmpEndIdx);
                        
                        % Linear interpolation across the peak.
                        replacementData         = interp1([neighbourFreqIndices(tmpStartIdx),neighbourFreqIndices(tmpEndIdx)],[tmpStartVal,tmpEndVal],indicesToReplace,'linear');
                        
                    case 'leaveSlope2'
                        replacementData               = justSlope(fittedStruct.b(chanIdx),fittedStruct.c(chanIdx),indicesToReplace);
                        
                end
                
                if strcmp(cfg.log,'yes')
                    replacementData = exp(replacementData);
                end
                
                % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
                fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), replacementData);

            end
        end
    case 'no'
        filtFreq        = zeros(1,2);

        % Check there are 2 or more peaks
        if length(peakFreq) < 2
            error('Only one peak found. Please set cfg.independentPeaks to yes');
        end
        
        % Find the bounds for the width
        peakLowerBound              = peakFreq - peakWidth;
        peakUpperBound              = peakFreq + peakWidth;
        widthFreqIndicesBound       = nearest(fftFreq,[min(peakLowerBound),max(peakUpperBound)]);
        
        % And for the specified neighbourhood.
        neighbourLowerBound         = peakFreq(1) - cfg.Neighwidth;
        neighbourUpperBound         = peakFreq(end) + cfg.Neighwidth;
        neighbourFreqIndicesBound   = nearest(fftFreq,[neighbourLowerBound,neighbourUpperBound]);
        neighbourFreqIndices        = neighbourFreqIndicesBound(1):neighbourFreqIndicesBound(end);
        
        if strcmp(cfg.log,'yes')
            % Get the chan avg data for those bounds.
            neighbourData               = log(scaledWelchPow(:,neighbourFreqIndices));
        elseif strcmp(cfg.log,'no')
            % Get the chan avg data for those bounds.
            neighbourData               = scaledWelchPow(:,neighbourFreqIndices);
        end
               
        % Guesses for fit (channel level)
        A           = zeros(size(peakFreq));
        x0          = zeros(size(peakFreq));
        g           = zeros(size(peakFreq));
        for peakIdx = 1:length(peakFreq)
            % Guesses for fit (peak level)
            % Amplitude
            % Floor
            A(peakIdx)              = peakProminence(peakIdx);
            
            % Centre
            x0(peakIdx)             = nearest(fftFreq,peakFreq(peakIdx));

            % Gamma. Sqrt to account for smoothing of estimate.
            g(peakIdx)              = nearest(fftFreq,peakWidth(peakIdx));
        end

        
        for chanIdx = 1:length(fftData(:,1))
            % Guesses for fit (channel level)
            neighbourChanData   = neighbourData(chanIdx,:);
            % Slope
            quarterLength   = round(length(neighbourChanData)/4);
            endVals         = mean(neighbourChanData((end - quarterLength):end));
            startVals       = mean(neighbourChanData(1:1+quarterLength));
            eighthLength    = round(quarterLength/2);
            b               = (endVals - startVals) / (neighbourFreqIndices(end - eighthLength) - neighbourFreqIndices(eighthLength));

            % Constant - What is c equal to for the slope go through the middle?
            middleValue     = mean([endVals,startVals]);
            c               = middleValue - (b * mean(x0));

            % Fit two peaks with slope to the channel data.
            fittedModel     = fit(neighbourFreqIndices', neighbourChanData', twoPeaksWithSlope,...
                                'StartPoint', [A(1), x0(1), g(1), A(2), x0(2), g(2), b, c],...
                                'Lower',[-inf, x0(1)-g(1), -inf, -inf, x0(2)-g(2), -inf, -inf, -inf],...
                                'Upper',[inf, x0(1)+g(1), inf, inf, x0(2)+g(2), inf, inf, inf]);

            fitCoeffValue   = coeffvalues(fittedModel);
            fitCoeffNames   = coeffnames(fittedModel);

            % Dynamic naming seems to be the best way...
            for coeffIdx = 1:length(fitCoeffNames)
                fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
            end
        end
        
        % Better estimate of x0
        x0(1)           = median(fittedStruct.x01);
        x0(2)           = median(fittedStruct.x02);
        
        
        % It might be worth using the maxG as a limit for the fitted
        % model. 
        medianG(1)        = median(abs(fittedStruct.g1));
        iqrG(1)            = iqr(abs(fittedStruct.g1));
        maxG(1)            = medianG(1) + 3*iqrG(1);
        minG(1)            = abs(medianG(1) - 3*iqrG(1));
        medianG(2)        = median(abs(fittedStruct.g2));
        iqrG(2)            = iqr(abs(fittedStruct.g2));
        maxG(2)            = medianG(2) + 3*iqrG(2);
        minG(2)            = abs(medianG(2) - 3*iqrG(2));
        
        for chanIdx = 1:length(fftData(:,1))
            % Guesses for fit (channel level)
            neighbourChanData   = neighbourData(chanIdx,:);
            % Slope
            quarterLength   = round(length(neighbourChanData)/4);
            endVals         = mean(neighbourChanData((end - quarterLength):end));
            startVals       = mean(neighbourChanData(1:1+quarterLength));
            eighthLength    = round(quarterLength/2);
            b               = (endVals - startVals) / (neighbourFreqIndices(end - eighthLength) - neighbourFreqIndices(eighthLength));

            % Constant - What is c equal to for the slope go through the middle?
            middleValue     = mean([endVals,startVals]);
            c               = middleValue - (b * mean(x0));

            % Fit two peaks with slope to the channel data.
            fittedModel     = fit(neighbourFreqIndices', neighbourChanData', twoPeaksWithSlope,...
                                'StartPoint', [A(1), x0(1), medianG(1), A(2), x0(2), medianG(2), b, c],...
                                'Lower',[-inf, x0(1), minG(1), -inf, x0(2), minG(2), -inf, -inf],...
                                'Upper',[inf, x0(1), maxG(1), inf, x0(2), maxG(2), inf, inf]);

            fitCoeffValue   = coeffvalues(fittedModel);
            fitCoeffNames   = coeffnames(fittedModel);

            % Dynamic naming seems to be the best way...
            for coeffIdx = 1:length(fitCoeffNames)
                fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
            end
        end
        
        maxG(1)            = max(abs(fittedStruct.g1));
        maxG(2)             = max(abs(fittedStruct.g2));
        maxG                = max(maxG);
        
        % Get the width of fitted lorentzians.
        x0Min                   = min([fittedStruct.x01(1);fittedStruct.x02(1)]);
        x0Max                   = max([fittedStruct.x01(1);fittedStruct.x02(1)]);
        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * maxG));
        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * maxG));

        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
        
        
        filtFreq(1,1:2)   = [fftFreq(indicesToReplace(1)),fftFreq(indicesToReplace(end))];

        if length(indicesToReplace) > length(neighbourData)
            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
        end
        
        for chanIdx = 1:length(fftData(:,1))
            switch cfg.method
                case 'removePeak'
                    % Get the lorentzian component.
                    peaksOnly             = justTwoPeaks(fittedStruct.A1(chanIdx),fittedStruct.x01(chanIdx),fittedStruct.g1(chanIdx),fittedStruct.A2(chanIdx),fittedStruct.x02(chanIdx),fittedStruct.g2(chanIdx),indicesToReplace);

                    % Remove it from the original data.
                    replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - peaksOnly;


                case 'leaveSlope'
                    % Get the peak on a slope
                    fittedModelData         = twoPeaksWithSlope(fittedStruct.A1(chanIdx),fittedStruct.x01(chanIdx),fittedStruct.g1(chanIdx),fittedStruct.A2(chanIdx),fittedStruct.x02(chanIdx),fittedStruct.g2(chanIdx),fittedStruct.b(chanIdx),fittedStruct.c(chanIdx),neighbourFreqIndices);
                    tmpStartVal             = fittedModelData(tmpStartIdx);
                    tmpEndVal               = fittedModelData(tmpEndIdx);

                    % Linear interpolation across the peak.
                    replacementData  = interp1([neighbourFreqIndices(tmpStartIdx),neighbourFreqIndices(tmpEndIdx)],[tmpStartVal,tmpEndVal],indicesToReplace,'linear');

            end
            
            if strcmp(cfg.log,'yes')
                replacementData = exp(replacementData);
            end
            
            % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
            fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), replacementData);
        end
end

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filteredTimeseries = ifft(fftData,[],2,'symmetric');

% Put it back into original structure.
filt            = data;
filt.trial{1}   = filteredTimeseries;

