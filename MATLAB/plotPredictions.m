%% Plotting predictions
% For a specific run
function plotPredictions(YPredNormalized, X, Y, V, downsampleFactor)
    minValueY = [0 0 0 0];
    maxValueY = [20/V 20/V 20/V 20/V];
    YPred = YPredNormalized .* (maxValueY - minValueY) + minValueY;
    YDownsampled = Y(15001:downsampleFactor:16000, :)*0.001;
    YPredDownsampled = YPred(15001:downsampleFactor:16000, :)*0.001;
    timeDownsampled = X(15001:downsampleFactor:16000, 1);
    
    figure;
    hold all
    
    for i = 1:4
        plot(timeDownsampled, YDownsampled(:, i), 'o', 'MarkerSize', 4);
        plot(timeDownsampled, YPredDownsampled(:, i), 'LineWidth', 2);
    end
    
    title ('Esterification reaction of methanol by acetic acid in a batch reactor');
    legend('C AA', 'C M', 'C MA', 'C W');
    xlabel('Time [s]');
    ylabel('Concentration [mol/L]');
    grid on;
    hold off;
end