%% Normalization of data
% Values used to normalize between 0.1 and 0.9 the inputs and outputs of the different neural networks.
function [XNormalized, YNormalized] = normalizeData(X, Y, V)
    % Range of values
    minValueX = [0 273 0 0 0 0 0];
    maxValueX = [1 373 15 20/V 20/V 20/V 20/V];
    minValueY = [0 0 0 0];
    maxValueY = [20/V 20/V 20/V 20/V];

    XNormalized = (X - minValueX) ./ (maxValueX - minValueX);

    YNormalized = (Y - minValueY) ./ (maxValueY - minValueY);
end
