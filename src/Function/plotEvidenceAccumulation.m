function plotEvidenceAccumulation(SelTrial, Tk, Ck, pp, ipp, ValueClasses, NumTrials, Threshold, LbClasses, CueClasses )

    cindex = Tk == SelTrial;

    [~, ClassIdx] = ismember(unique(Ck(cindex)), CueClasses);

    GreyColor = [150 150 150]/255;
    LineColors = {'b', 'g', 'r'};
    
    hold on;
    
    % Plotting raw probabilities
    plot(pp(cindex, 1), 'o', 'Color', GreyColor);

    % Plotting accumulutated evidence
    plot(ipp(cindex), 'k', 'LineWidth', 2);

    % Plotting actual target class
    yline(ValueClasses(ClassIdx), LineColors{ClassIdx}, 'LineWidth', 5);

    % Plotting 0.5 line
    yline(0.5, '--k');

    % Plotting thresholds
    yline(Threshold, 'k', 'Th_{1}');
    yline(1-Threshold, 'k', 'Th_{2}');
    hold off;

    grid on;
    ylim([0 1]);
    xlim([1 sum(cindex)]);
    legend('raw prob', 'integrated prob');
    ylabel('probability/control')
    xlabel('sample');
    title(['Trial ' num2str(SelTrial) '/' num2str(NumTrials) ' - Class ' LbClasses ' (' num2str(CueClasses(ClassIdx)) ')']);

end