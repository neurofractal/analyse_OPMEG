function [filt] = ft_preproc_dft_remove_gauss(cfg, data)
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
%   cfg.method              = 'leaveSlopeMinusGauss', 'removeGauss' 
%                               'removeLorentzian', or 
%                               'leaveSlopeMinusLorentzian'
%   cfg.independentPeaks    = 'yes' or 'no'. Whether to model peaks 
%                               independently or together.
%   cfg.log                 = 'yes' or 'no'. Whether to log the PSD before
%                               fitting.

%% Extract info from input structure
timeSeries          = data.trial{1};
nsamples            = length(timeSeries);
samplingFreq        = data.fsample;
fftFreq             = samplingFreq * linspace(0, 1, nsamples);

% Run the fft
fftData             = fft(timeSeries,nsamples,2);
avgData             = median(abs(fftData),1);

%% Identify peaks in power in PSD
% Calculate PSD in smaller epochs. 
cfgPSD              = [];
cfgPSD.channel      = 'all';
cfgPSD.trial_length = 3;
cfgPSD.method       = 'tim'; % Breaks into epochs
cfgPSD.foi          = cfg.foi;
cfgPSD.plot         = 'no';
[pow,freq,~]        = ft_opm_psd(cfgPSD,data);

% Use those to produce a smoother PSD output for peak identification.
powMed              = median(pow(:,:,:),3);
powMed              = median(pow(:,:),2);
powMedLog           = log(powMed);

% And interpret to the same index as the main fft (makes things easier)
if strcmp(cfg.log,'yes')
    smoothPow           = interp1(freq,powMedLog,fftFreq);
elseif strcmp(cfg.log,'no')
    smoothPow           = interp1(freq,powMed,fftFreq);
end

% Tidy up
clear pow powMed powMedLog cfgPSD freq samplingFreq cfgPSD timeSeries...
nsamples

% Select pow data for the foi.
foiIdx            	= ~(fftFreq < cfg.foi(1) | fftFreq > cfg.foi(2));
smoothPowFOI        = smoothPow(foiIdx);
freqFOI             = fftFreq(foiIdx);

% The user will be asked if they are satisfied with peak identification.
noPeaks             = true;
firstRun            = true;
while noPeaks || changePeaks
    % Get a staring point for peak prominence in the first run.
    if firstRun
        [~,~,~,initialPeakPromList]   = findpeaks(smoothPowFOI,freqFOI);
        minPeakProminence           = max(initialPeakPromList);
    end
    % Identify peaks based on prominence. Extract halfheight width.
    if strcmp(cfg.independentPeaks,'yes')
        [peakVal,peakFreq,peakWidth,peakProminence]  = findpeaks(smoothPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'MinPeakDistance',cfg.minPeakDistance,'WidthReference','halfprom');
    elseif strcmp(cfg.independentPeaks,'no')
        [peakVal,peakFreq,peakWidth,peakProminence]  = findpeaks(smoothPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'WidthReference','halfprom');
    end

    if isempty(peakVal)
        minPeakProminence = minPeakProminence / 2;
    else
        noPeaks = false;
        if strcmp(cfg.independentPeaks,'yes')
            findpeaks(smoothPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'Annotate','extents','MinPeakDistance',cfg.minPeakDistance,'WidthReference','halfprom');
        elseif strcmp(cfg.independentPeaks,'no')
            findpeaks(smoothPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'Annotate','extents','WidthReference','halfprom');
        end
        
    end

    if ~noPeaks
        % Ask if the user wants to change the number of peaks
        diffPeaks           = 3;
        while ~any(eq(diffPeaks,[0 1 2]))
            diffPeaks = input('To continue, enter 0\nFor less peaks, enter 1\nFor more peaks, enter 2:\n');
        end

        % There is probably a better way, but this works for now.
        if diffPeaks == 1
            if minPeakProminence == max(initialPeakPromList)
                disp('Already at minimum peaks. Enter either 2 or 0');
                close(gcf);
                changePeaks     = true;
            else
                minPeakProminence = minPeakProminence * 1.25;
                close(gcf);
                changePeaks     = true;
            end
        elseif diffPeaks == 2
            minPeakProminence = minPeakProminence / 1.25;
            close(gcf);
            changePeaks     = true;
        else
            changePeaks = false;
            close(gcf);
        end
    end
    if firstRun
        firstRun = false;
    end
end

%% Remove Gaussian component from amplitude in spectrum
switch cfg.independentPeaks
    case 'yes'
        for peakIdx = 1:length(peakFreq)
            % Find the bounds for the width
            peakLowerBound              = peakFreq(peakIdx) - peakWidth(peakIdx);
            peakUpperBound              = peakFreq(peakIdx) + peakWidth(peakIdx);
            widthFreqIndicesBound       = nearest(fftFreq,[peakLowerBound,peakUpperBound]);

            % And for the specified neighbourhood.
            neighbourLowerBound         = peakFreq(peakIdx) - cfg.Neighwidth;
            neighbourUpperBound         = peakFreq(peakIdx) + cfg.Neighwidth;
            neighbourFreqIndicesBound   = nearest(fftFreq,[neighbourLowerBound,neighbourUpperBound]);
            neighbourFreqIndices        = neighbourFreqIndicesBound(1):neighbourFreqIndicesBound(end);
            
            % Find the bounds for the specified neighbourhood.
            indivNeighbourLowerBound         = peakFreq(peakIdx) - cfg.Neighwidth;
            indivNeighbourUpperBound         = peakFreq(peakIdx) + cfg.Neighwidth;
            indivNeighbourFreqIndices        = nearest(fftFreq,[indivNeighbourLowerBound,indivNeighbourUpperBound]);
            
            if strcmp(cfg.log,'yes')
                % Get the chan avg data for those bounds.
                neighbourData               = log(abs(avgData(neighbourFreqIndices)));
                % Get the chan avg data for those bounds.
                indivNeighbourData               = log(abs(avgData(neighbourFreqIndices)));
            elseif strcmp(cfg.log,'no')
                % Get the chan avg data for those bounds.
                neighbourData               = abs(avgData(neighbourFreqIndices));
                % Get the chan avg data for those bounds.
                indivNeighbourData               = abs(avgData(neighbourFreqIndices));
            end

            % Guesses for fit (peak level)
            % Amplitude
            A               = peakProminence(peakIdx);
            
            % Centre
            x0          = nearest(fftFreq,peakFreq(peakIdx));

            % Gamma. Sqrt to account for smoothing of estimate.
            g                               = sqrt(nearest(fftFreq,peakWidth(peakIdx)));
            
            % Guesses for fit (channel level)
            % Slope
            endVals             = mean(neighbourData((end - (round(g/4))):end));
            startVals           = mean(neighbourData(1:1+(round(g/4))));
            b                   = (endVals - startVals) / (neighbourFreqIndicesBound(end) - neighbourFreqIndicesBound(1));

            % Constant
            % What is c equal to for the slope go through the middle?
            middleValue     = mean([endVals,startVals]);
            c               = middleValue - (b * median(neighbourFreqIndicesBound));
            

            % Define function for Gaussian/Lorentzian with and without
            % slopes.
            gaussWithSlope = @(A, x0, s, b, c, x)...
                (A*exp(-2*((x-x0)/s).^2)) + (b*x) + c;

            gauss = @(A, x0, s, x)...
                (A*exp(-2*((x-x0)/s).^2));

            slope = @(b, c, x)...
                (b*x) + c;

            lorentzianWithSlope = @(A, x0, g, b, c, x)...
                ((A .* (g.^2))./((g.^2) + ((x - x0).^2))) + (b.*x) + c;

            lorentzian = @(A, x0, g, x)...
                ((A .* (g.^2))./((g.^2) + ((x - x0).^2)));
            
            switch cfg.method
                case 'leaveSlopeMinusGauss'
                    % Fit a gauss with slope to the data.
                    fittedAvgModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlope,...
                                                'StartPoint', [A, x0, g, b, c]);

                    bestGuess           = gaussWithSlope(A, x0, g, b, neighbourFreqIndices');
                    
                case 'removeGauss'
                    % Fit a gauss with slope to the data.
                    fittedAvgModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlope,...
                                                'StartPoint', [A, x0, g, b, c]);
                    
                    bestGuess           = gaussWithSlope(A, x0, g, b, neighbourFreqIndices');
                    
                case 'removeLorentzian'
                    % Fit a Lorentzian with slope
                    fittedAvgModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlope,...
                        'StartPoint', [A, x0, g, b, c]);
                    
                    bestGuess = lorentzianWithSlope(A, x0, g, b, c, neighbourFreqIndices');
                    
                case 'leaveSlopeMinusLorentzian'
                    % Fit a lortenzian with slope to the data.
                    fittedAvgModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlope,...
                        'StartPoint', [A, x0, g, b, c]);
                    
                    fftFreqIndices = 1:length(fftFreq);
                    bestGuess = lorentzianWithSlope(A, x0, g, b, c, fftFreqIndices');
            end
            
% %             See how good the fit was.
%             slopeOnly = slope(b, c, fftFreqIndices');
%             hold on
%             plot(neighbourFreqIndices, neighbourData);
%             plot(fftFreqIndices,bestGuess);
%             plot(fittedAvgModel);
%             plot(slopeOnly);
%             plot(fftFreqIndices,smoothPow + 10);
%             hold off
%             
            % Redefine functions with constants for peak and width.
            gaussWithSlopeC = @(A, b, c, x)...
                (A*exp(-2*((x-fittedAvgModel.x0)/fittedAvgModel.s).^2)) + (b*x) + c;

            gaussC = @(A, x)...
                (A*exp(-2*((x-fittedAvgModel.x0)/fittedAvgModel.s).^2));

            lorentzianWithSlopeC = @(A, b, c, x)...
                ((A * (fittedAvgModel.g^2))./((fittedAvgModel.g^2) + ((x - fittedAvgModel.x0).^2))) + (b*x) + c;

            lorentzianC = @(A, x)...
                ((A * (fittedAvgModel.g^2))./((fittedAvgModel.g.^2) + ((x - fittedAvgModel.x0).^2)));
            
            
            % Fit channels independently
            for chanIdx = 1:length(fftData(:,1))
                if strcmp(cfg.log,'yes')
                    % Get the data being fit.
                    peakWidthData = log(abs(fftData(chanIdx,widthFreqIndicesBound(1):widthFreqIndicesBound(end))));
                    neighbourData   = log(abs(fftData(chanIdx,neighbourFreqIndices)));
                elseif strcmp(cfg.log,'no')
                    % Get the data being fit.
                    peakWidthData = abs(fftData(chanIdx,widthFreqIndicesBound(1):widthFreqIndicesBound(end)));
                    neighbourData   = abs(fftData(chanIdx,neighbourFreqIndices));
                end
                
                switch cfg.method
                    case 'leaveSlopeMinusGauss'
                        % Fit a gauss with slope to the data.
                        fittedModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlopeC,...
                                                    'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of 3 std.dev of the fitted Gauss
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (3 * fittedAvgModel.s));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (3 * fittedAvgModel.s));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                        % Get just the slope
                        replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);

                    case 'removeGauss'
                        % Fit a gauss with slope to the data.
                        fittedModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlopeC,...
                                                    'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of 3 std.dev of the fitted Gauss
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (3 * fittedAvgModel.s));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (3 * fittedAvgModel.s));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                        % Find the Gauss component.
                        onlyGauss               = gaussC(fittedModel.A, indicesToReplace);

                        % Remove the Gauss from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyGauss;

                    case 'removeLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlopeC,...
                            'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of 6 gamma of fitted lorentzian.
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (6 * abs(fittedAvgModel.g)));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (6 * abs(fittedAvgModel.g)));

                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end

                        % Get the lorentzian component.
                        onlyLorentz             = lorentzianC(fittedModel.A, indicesToReplace);

                        % Remove it from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;

                    case 'leaveSlopeMinusLorentzian'
                        % Fit a lortenzian with slope to the data.
                        fittedModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlopeC,...
                            'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of 6 gamma of the fitted Lorentz
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (6 * abs(fittedAvgModel.g)));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (6 * abs(fittedAvgModel.g)));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                        bestGuess       = lorentzianWithSlopeC(fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c, neighbourFreqIndices');

                        % Get just the slope
                        roughSlope                      = slope(fittedModel.b,fittedModel.c, indicesToReplace);
                        neighbourDataWithReplacement    = neighbourData;
                        replacementIndices              = ismember(neighbourFreqIndices, indicesToReplace);
                        neighbourDataWithReplacement(replacementIndices) = roughSlope;
                        fittedSlope                     = fit(neighbourFreqIndices',neighbourDataWithReplacement',slope,'StartPoint',[fittedAvgModel.b, fittedAvgModel.c]);
                        replacementData                 = fittedSlope(indicesToReplace);
                        replacementData                 = replacementData';
                        % The slope will always be slightly too low.
%                         hold on
%                         plot(neighbourFreqIndices,neighbourData)
%                         plot(indicesToReplace,replacementData);

%                         
%                         plot(fittedSlope,neighbourFreqIndices',tmp')
%                         plot(indicesToReplace,replacementData2);
                        
                end
                
%                 % debug plots
%                 subplot(3,1,1);
%                 hold on
%                 plot(neighbourFreqIndices, neighbourData)
%                 plot(fittedModel)
%                 hold off
%                 subplot(3,1,2);
%                 hold on
%                 plot(indicesToReplace,replacementData)
%                 plot(fittedModel)
%                 hold off
%                 subplot(3,1,3);
%                 hold on
%                 plot(fittedModel)
%                 plot(bestGuess)
%                 hold off
                
                % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
                if strcmp(cfg.log,'yes')
                    fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), exp(replacementData));
                elseif strcmp(cfg.log,'no')
                    fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), replacementData);
                end
                
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
            peakWidthData               = log(abs(avgData(widthFreqIndicesBound(1):widthFreqIndicesBound(end))));
            neighbourData               = log(abs(avgData(neighbourFreqIndices)));
        elseif strcmp(cfg.log,'no')
            % Get the chan avg data for those bounds.
            peakWidthData               = abs(avgData(widthFreqIndicesBound(1):widthFreqIndicesBound(end)));
            neighbourData               = abs(avgData(neighbourFreqIndices));
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
        slope = @(b, c, x)...
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
            tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax)); % Magic number
            tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax)); % Magic number
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
            tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
            tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));
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
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));

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
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));

                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end
                        
                        % Get just the slope
                        replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);
                        
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
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));
                        
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
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));                        
                        
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end
                        
                        % Get just the slope
                        replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);
                        
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

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filteredTimeseries = ifft(fftData,[],2,'symmetric');

% Put it back into original structure.
filt            = data;
filt.trial{1}   = filteredTimeseries;

