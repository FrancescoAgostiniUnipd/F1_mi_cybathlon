%% Fisher Scor of a set of runs
%   PSD,
%   Runs, Different classes inside PSD to see
%   Rk  , Sample appartain to Run #
%   Ck  , Cue is classified 771, 773, 0
%   classId, 2 classes to see index
function score = FisherScore(PSD,Runs, Rk,Ck,classId)
    % disp(classId);
    NumFreqs = size(PSD, 2);
    NumChans = size(PSD, 3);
    NumClasses = length( classId );
    NumRuns = length(Runs);
    score = nan(NumFreqs, NumChans, NumRuns);
    % FS2 = nan(NumFreqs*NumChans, NumRuns);
    for rId = 1:NumRuns
        rindex = Rk == Runs(rId); 
        
        cmu    = nan(NumFreqs, NumChans, 2);
        csigma = nan(NumFreqs, NumChans, 2);
        
        for cId = 1:NumClasses
            cindex = rindex & Ck == classId(cId);
            cmu(:, :, cId) = squeeze(mean(PSD(cindex, :, :)));
            csigma(:, :, cId) = squeeze(std(PSD(cindex, :, :)));
        end
        
        score(:, :, rId) = abs(cmu(:, :, 2) - cmu(:, :, 1)) ./ sqrt( ( csigma(:, :, 1).^2 + csigma(:, :, 2).^2 ) );
    end
end