%% plot Fisher Score
%   fisherScore, score  of fisher to plot
%   NumRun      , # of score to plot
%   freqs       , array of select frequency
%   channelLb   , array of Channels
function  plotFisherScore(fisherScore, NumRuns,freqs,channelLb)
    imagesc(fisherScore(:, :, NumRuns)');
    axis square;
    NumFreqs = size(freqs);
    NumChans = size(channelLb,2);
    set(gca, 'XTick', 1:NumFreqs);
    set(gca, 'XTickLabel', freqs);
    set(gca, 'YTick', 1:NumChans);
    set(gca, 'YTickLabel', channelLb);
    xtickangle(-90);
    % colorbar;
    title(['Calibration run ' num2str(NumRuns)]);
end