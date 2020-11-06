function [filt, filteredFrequencies] = ft_preproc_dft_remove_gauss(cfg, data)
% Function to find peaks in spectrum, model their shape and remove them 
% before transforming back into time domain. 
% 
% Function is based on Search Results ft_preproc_dftfilter. In particular,
% the spectrum interpolation method (Leske & Dalal, 2019, NeuroImage 189,
% doi: 10.1016/j.neuroimage.2019.01.026).
% 
% Overview of function
% 1)   Data transformed into the frequency domain via a discrete Fourier 
%       transform (DFT). 
% 2)   Peaks in the specified frequency range are identified and some
%       characteristics extracted about them.
% 3)   The characteristics are used to model them onto data from each
%       channel.
% 4)   A specified method is used to remove the peaks.
% 4)   The signal is transformed back into the time domain via inverse DFT
%       (iDFT). 
%
% Use as
%   [filt] = ft_preproc_dft_remove_gause(cfg, data)
%   where
%   data                    = Fieldtrip data structure with continuous
%                               data. Usually raw data.
%   cfg                     = Structure containing the following:
%   cfg.foi                 = [start end], freq range to remove peaks 
%                               between e.g. [40 60];
%   cfg.Neighwidth          = width (in Hz) of peaks to evaluate, e.g. 2;
%   cfg.minPeakDistance     = minimum distance in Hz between peaks.
%   cfg.method              = 'leaveSlope' or 'removePeak' 
% 	cfg.peakShape           = 'Guassian' or 'Lorentzian'
%   cfg.independentPeaks    = 'yes' or 'no'. Whether to model peaks 
%                               independently or together.
%   cfg.log                 = 'yes' or 'no'. Whether to log the PSD before
%                               fitting.
%   cfg.strength            = multiplier of g, the variable that determines 
%                             distribution width.

%% Extract info from input structure
timeSeries          = data.trial{1};
nsamples            = length(timeSeries);
samplingFreq        = data.fsample;
fftFreq             = samplingFreq * linspace(0, 1, nsamples);

% Run the fft
fftData             = fft(timeSeries,nsamples,2);

% Get real pow
% fftPow              = abs(fftData);
% fftPow              = fftPow(:,1:nsamples/2+1);
fftFreq             = fftFreq(1:nsamples/2+1);
% avgFftPow           = median(fftPow,1);
welchPow            = [];
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
% scaledWelchPow  = scaledWelchPow/2;
% scaledWelchPow  = scaledWelchPow*(window'*window);
% scaledWelchPow  = sqrt(scaledWelchPow);  
avgWelchPow     = median(scaledWelchPow,1);

% % Interp welch to fft freq resolution
% interpWelchPow  = interp1(welchFreq,avgWelchPow,fftFreq);

% Debug plot
% plot(fftFreq,avgFftPow,fftFreq,avgWelchPow,'r');
% set(gca, 'YScale', 'log')

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
    case 'Guassian'
        peakWithSlope   = @(A, x0, g, b, c, x)...
                            (A*exp(-2*((x-x0)/g).^2)) + (b*x) + c;
        justPeak        = @(A, x0, g, x)...
                            (A*exp(-2*((x-x0)/g).^2));
    case 'Lorentzian'
        peakWithSlope   = @(A, x0, g, b, c, x)...
                            ((A .* (g.^2))./((g.^2) + ((x - x0).^2))) + (b.*x) + c;
        justPeak        = @(A, x0, g, x)...
                            ((A .* (g.^2))./((g.^2) + ((x - x0).^2)));
end

%% Remove Gaussian component from amplitude in spectrum
switch cfg.independentPeaks
    case 'yes'
        filteredFrequencies        = zeros(length(peakFreq),2);
        for peakIdx = 1:length(peakFreq)
            %% Find the neighbourhood data around the peak
            neighbourLowerBound         = peakFreq(peakIdx) - (cfg.Neighwidth);
            neighbourUpperBound         = peakFreq(peakIdx) + (cfg.Neighwidth);
            neighbourFreqIndicesBound   = nearest(fftFreq,[neighbourLowerBound,neighbourUpperBound]);
            neighbourFreqIndices        = neighbourFreqIndicesBound(1):neighbourFreqIndicesBound(end);
            
%             % And the frequencies
%             neighbourFreq               = fftFreq(neighbourFreqIndices);
            
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
                                    'Lower',[-inf, x0, -inf, -inf, -inf],...
                                    'Upper',[inf, x0, inf, inf, inf]);

                fitCoeffValue   = coeffvalues(fittedModel);
                fitCoeffNames   = coeffnames(fittedModel);
                
                % Dynamic naming seems to be the best way...
                for coeffIdx = 1:length(fitCoeffNames)
                    fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
                end
            end
            
            % A central measure may work better here...
            maxG        = median(abs(fittedStruct.g));
            % Get the width of cfg.strength gamma of fitted lorentzian.
            tmpStartIdx             = nearest(neighbourFreqIndices,(x0 - (cfg.strength * maxG)));
            tmpEndIdx               = nearest(neighbourFreqIndices,(x0 + (cfg.strength * maxG)));

            indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
            
            filteredFrequencies(peakIdx,1:2)   = [fftFreq(indicesToReplace(1)),fftFreq(indicesToReplace(end))];
            if length(indicesToReplace) > length(neighbourData(chanIdx,:))
                ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
            end
            
            for chanIdx = 1:length(fftData(:,1))
                switch cfg.method
                    case 'removePeak'
                        % Get just the peak
                        peakOnly                    = justPeak(fittedStruct.A(chanIdx),fittedStruct.x0(chanIdx),fittedStruct.g(chanIdx),indicesToReplace);

                        % Remove it from the original data.
                        replacementData  = neighbourData(chanIdx,tmpStartIdx:tmpEndIdx) - peakOnly;

                    case 'leaveSlope'
                        % Get the peak on a slope
                        fittedModelData         = peakWithSlope(fittedStruct.A(chanIdx),fittedStruct.x0(chanIdx),fittedStruct.g(chanIdx),fittedStruct.b(chanIdx),fittedStruct.c(chanIdx),neighbourFreqIndices);
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
    case 'no'
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
            neighbourData               = log(scaledWelchPow(neighbourFreqIndices));
        elseif strcmp(cfg.log,'no')
            % Get the chan avg data for those bounds.
            neighbourData               = scaledWelchPow(neighbourFreqIndices);
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
            g(peakIdx)              = sqrt(nearest(fftFreq,peakWidth(peakIdx)));
        end
        
        % Slope
        endVals             = mean(neighbourData((end - (round(g/4))):end));
        startVals           = mean(neighbourData(1:1+(round(g/4))));
        b                   = (endVals - startVals) / (neighbourFreqIndicesBound(end) - neighbourFreqIndicesBound(1));

        % Constant
        middleValue     = mean([endVals,startVals]);
        c               = middleValue - (b * median(neighbourFreqIndicesBound));
        
        % Define slope
        justSlope = @(b, c, x)...
                (b*x) + c;

        % Fit to the avg of channel data to get best guess.
        if length(peakFreq) == 2
            twoLorentzianWithSlope = @(A1, x01, g1, A2, x02, g2, b, c, x)...
                    ((A1 .* (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2))) + (b.*x) + c;

            twoLorentzian = @(A1, x01, g1, A2, x02, g2, x)...
                    ((A1 * (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2)));

            fittedAvgModel     = fit(neighbourFreqIndices', neighbourData', twoLorentzianWithSlope,...
                                'StartPoint', [A(1), x0(1), g(1), A(2), x0(2), g(2), b, c]);

            % Get the width of fitted lorentzians.
            x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02]);
            x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02]);
            gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2]));
            tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax)); % Magic number
            tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax)); % Magic number
            indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

            if length(indicesToReplace) > length(neighbourData)
                ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
            end

        elseif length(peakFreq) == 3
            threeLorentzianWithSlope = @(A1, x01, g1, A2, x02, g2, A3, x03, g3, b, c, x)...
                    ((A1 .* (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2))) + ((A3 * (g3.^2))./((g3.^2) + ((x - x03).^2))) + (b.*x) + c;

            threeLorentzian = @(A1, x01, g1, A2, x02, g2, A3, x03, g3, x)...
                    ((A1 .* (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2))) + ((A3 * (g3.^2))./((g3.^2) + ((x - x03).^2)));

            fittedAvgModel     = fit(neighbourFreqIndices', neighbourData', threeLorentzianWithSlope,...
                                'StartPoint', [A(1), x0(1), g(1), A(2), x0(2), g(2), A(3), x0(3), g(3), b, c]);
            
            % Get the width of fitted lorentzians.
            x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
            x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
            gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2,fittedAvgModel.g3]));
            tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
            tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));
        end

        % Now apply per channel with some constants from the avg.  
        for chanIdx = 1:length(fftData(:,1))
            if strcmp(cfg.log,'yes')
                peakWidthData = log(abs(fftData(chanIdx,widthFreqIndicesBound(1):widthFreqIndicesBound(end))));
                neighbourData   = log(abs(fftData(chanIdx,neighbourFreqIndices)));
            elseif strcmp(cfg.log,'no')
                peakWidthData = abs(fftData(chanIdx,widthFreqIndicesBound(1):widthFreqIndicesBound(end)));
                neighbourData   = abs(fftData(chanIdx,neighbourFreqIndices));
            end
            

            if length(peakFreq) == 2
                 twoLorentzianWithSlopeC = @(A1, A2, b, c, x)...
                    ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + (b*x) + c;

                twoLorentzianC = @(A1, A2, x)...
                    ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2)));

                switch cfg.method
                    case 'leaveSlopeMinusGauss'
                        
                    case 'removeGauss'
                        
                    case 'removeLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel     = fit(neighbourFreqIndices', neighbourData', twoLorentzianWithSlopeC,...
                                            'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of fitted lorentzians.
                        x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02]);
                        x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02]);
                        gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2]));
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));

                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end

                        % Get the lorentzian component.
                        onlyLorentz             = twoLorentzianC(fittedModel.A1, fittedModel.A2, indicesToReplace);

                        % Remove it from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;

                        bestGuess       = twoLorentzianWithSlopeC(fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c,neighbourFreqIndices');

                    case 'leaveSlopeMinusLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel     = fit(neighbourFreqIndices', neighbourData', twoLorentzianWithSlopeC,...
                                            'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of fitted lorentzians.
                        x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02]);
                        x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02]);
                        gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2]));
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));

                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end
                        
                        % Get just the slope
                        replacementData         = justSlope(fittedModel.b,fittedModel.c, indicesToReplace);
                        
                        bestGuess       = twoLorentzianWithSlopeC(fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c,neighbourFreqIndices');
                end
                
            elseif length(peakFreq) == 3
                threeLorentzianWithSlopeC = @(A1, A2, A3, b, c, x)...
                        ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + ((A3 * (fittedAvgModel.g3^2))./((fittedAvgModel.g3^2) + ((x - fittedAvgModel.x03).^2))) + (b*x) + c;
                
                threeLorentzianC = @(A1, A2, A3, x)...
                        ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + ((A3 * (fittedAvgModel.g3^2))./((fittedAvgModel.g3^2) + ((x - fittedAvgModel.x03).^2)));

                
                
                switch cfg.method
                    case 'leaveSlopeMinusGauss'
                        
                    case 'removeGauss'
                        
                    case 'removeLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel     = fit(neighbourFreqIndices', neighbourData', threeLorentzianWithSlopeC,...
                                    'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of fitted lorentzians.
                        x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                        x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                        gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2,fittedAvgModel.g3]));
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));
                        
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end

                        % Get the lorentzian component.
                        onlyLorentz             = threeLorentzianC(fittedModel.A1, fittedModel.A2, fittedAvgModel.A3, indicesToReplace);

                        % Remove it from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;

                        bestGuess       = threeLorentzianWithSlopeC(fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c, neighbourFreqIndices');

                    case 'leaveSlopeMinusLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel     = fit(neighbourFreqIndices', neighbourData', threeLorentzianWithSlopeC,...
                                    'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of fitted lorentzians.
                        x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                        x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                        gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2,fittedAvgModel.g3]));
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));                        
                        
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end
                        
                        % Get just the slope
                        replacementData         = justSlope(fittedModel.b,fittedModel.c, indicesToReplace);
                        
                        bestGuess       = threeLorentzianWithSlopeC(fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c, neighbourFreqIndices');
                end
            end
%             % debug plots
%             subplot(3,1,1);
%             hold on
%             plot(neighbourFreqIndices, neighbourData)
%             plot(fittedModel)
%             hold off
%             subplot(3,1,2);
%             hold on
%             plot(indicesToReplace,replacementData)
%             plot(fittedModel)
%             hold off
%             subplot(3,1,3);
%             hold on
%             plot(fittedModel)
%             plot(bestGuess)
%             hold off

            % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
            if strcmp(cfg.log,'yes')
                fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), exp(replacementData));
            elseif strcmp(cfg.log,'no')
                fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), replacementData);
            end
        end
end

%% Plot the result
% % Original Data
% neighbourFreq       = fftFreq(neighbourFreqIndices);
% 
% neighbourFftPow     = welchPow(:,neighbourFreqIndices);
% errorOrigFftPow     = std(neighbourFftPow,[],1)./sqrt(size(neighbourFftPow,1));
% CI95OrigFftPow      = bsxfun(@plus, avgWelchPow(neighbourFreqIndices)', bsxfun(@times, [-1  1]*1.96, errorOrigFftPow'))';   % 95% Confidence Intervals
% 
% % Replaced Data
% newFftPow           = abs(fftData(:,indicesToReplace));
% avgNewFftPot        = median(newFftPow,1);
% replacedFreq        = fftFreq(indicesToReplace);
% errorNewFftPow     = std(newFftPow,[],1)./sqrt(size(newFftPow,1));
% CI95NewFftPow      = bsxfun(@plus, avgNewFftPot', bsxfun(@times, [-1  1]*1.96, errorNewFftPow'))';   % 95% Confidence Intervals
% 
% % Difference between them
% replacedOrigFftPow   = fftPow(:,indicesToReplace);
% differencesFftPow   = newFftPow - replacedOrigFftPow ;
% avgDifferencesFftPow    = median(differencesFftPow,1);
% errorDifferencesFftPow     = std(differencesFftPow,[],1)./sqrt(size(differencesFftPow,1));
% CI95DifferencesFftPow      = bsxfun(@plus, avgDifferencesFftPow', bsxfun(@times, [-1  1]*1.96, errorDifferencesFftPow'))';   % 95% Confidence Intervals
% 
% 
% % Plot result.
% figure
% hold on
% 
% % Original Data
% fill([neighbourFreq';flipud(neighbourFreq')]',[CI95OrigFftPow(1,:)';CI95OrigFftPow(2,:)']',[.9 .9 .9],'linestyle','none');
% line(neighbourFreq,avgFftPow(neighbourFreqIndices))
% 
% % New Data
% fill([replacedFreq';flipud(replacedFreq')]',[CI95NewFftPow(1,:)';CI95NewFftPow(2,:)']',[.9 .9 .9],'linestyle','none');
% line(replacedFreq,avgNewFftPot)
% 
% % Difference data
% fill([replacedFreq';flipud(replacedFreq')]',[CI95DifferencesFftPow(1,:)';CI95DifferencesFftPow(2,:)']',[.9 .9 .9],'linestyle','none');
% line(replacedFreq,avgDifferencesFftPow)
% 
% set(gca, 'YScale', 'log')
% 

% % Subtracted Data
% fill([replacedFreq';flipud(replacedFreq')]',[CI95SubtractedData(1,:)';CI95SubtractedData(2,:)']',[.9 .9 .9],'linestyle','none');
% line(replacedFreq,avgSubtractedData)
% 
% 

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filteredTimeseries = ifft(fftData,[],2,'symmetric');

% Put it back into original structure.
filt            = data;
filt.trial{1}   = filteredTimeseries;

