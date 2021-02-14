%% function plotEvidenceAccumulation(SelTrial, Tk, Ck, pp, ipp, ValueClasses, NumTrials, Threshold, LbClasses, CueClasses )
%   SelTrial,   Select Trial to visualize
%   Tk,         Trial logic vector
%   Ck,         Cue
%   pp,         Confidence parameter
%   ipp,        is parameter of integrated of pp: ipp(sId) = prev_ipp.*alpha + curr_pp.*(1-alpha)
%   ValueClasses,
%   NumTrials,  
%   Threshold,  
%   LbClasses,  
%   CueClasses  
function plotEvidenceAccumulation(SelTrial, Tk, Ck, pp, ipp, CueType, NumTrials, Threshold )

    cindex = Tk == SelTrial;
    CueClasses    = [771 783 773];

    [~, ClassIdx] = ismember(CueType(SelTrial), CueClasses);

    if ClassIdx == 0
        ClassIdx =783;
    end

    GreyColor = [150 150 150]/255;
    LineColors = {'b', 'g', 'r'};
    LbClasses     = {'both feet', 'rest', 'both hands'};
    ValueClasses  = [1 0.5 0];
    
    hold on;
    
    % Plotting raw probabilities
    plot(pp(cindex, 1), 'o', 'Color', GreyColor);

    % Plotting accumulutated evidence
    plot(ipp(cindex), 'k', 'LineWidth', 2);
    
    % Plotting actual target class
    yline(ValueClasses(ClassIdx), LineColors{ClassIdx}, 'LineWidth', 5);
    
    % Plotting thresholds
    yline(Threshold, 'k', 'Th_{1}');
    yline(1-Threshold, 'k', 'Th_{2}');

    % Plotting 0.5 line
    yline(0.5, '--k');
    hold off;

    grid on;
    ylim([0 1]);
    xlim([1 sum(cindex)]);
    legend('raw prob', 'integrated prob', 'Actual Target','Thr 1','Thr 2');
    ylabel('probability/control')
    xlabel('sample');
    title(['Trial ' num2str(SelTrial) '/' num2str(NumTrials) ' - Class ' LbClasses{ClassIdx} ' (' num2str(CueType(SelTrial)) ')']);

end