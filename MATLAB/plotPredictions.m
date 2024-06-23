%% Plotting predictions
% For a specific run
function plotPredictions(YPredNormalized, X, Y, run, V, downsampleFactor)
    if run == 1
        runIndex = 1;
    else
        runIndex = (run - 1) * 1000 + 1;
    end

    endIndex = run * 1000;

    minValueY = [0 0 0 0];
    maxValueY = [20/V 20/V 20/V 20/V];
    YPred = YPredNormalized .* (maxValueY - minValueY) + minValueY;
    YDownsampled = Y(runIndex:downsampleFactor:endIndex, :)*0.001;
    YPredDownsampled = YPred(runIndex:downsampleFactor:endIndex, :)*0.001;
    timeDownsampled = X(runIndex:downsampleFactor:endIndex, 1);
    
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