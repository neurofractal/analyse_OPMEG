function [yMean, yCI95] = get_confidence(y)
%
% Function to calculate mean and 95% confidence intervals for a matrix of
% data

% Input:    - y (organised as matrix with size N*M, where N is dimension
%           for averaging)
%
% Outputs:  - yMean = the mean of the input
%           - yCI95 = upper and lower limits for 95% confidence intervals

N = size(y,1);                                      
yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ?x?

ySEM = std(y)/sqrt(N);                 % Compute ?Standard Error Of The Mean? Of All Experiments At Each Value Of ?x?
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ?x?
end