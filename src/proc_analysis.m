%%Proc Script

clearvars;
% Adding functions class
addpath('./Function/');
addpath('./Util/eeglab/sigprocfunc/');
addpath('./Util/eeglab/adminfunc/');
addpath('./Util/eeglab/guifunc/');
%% Load data

% crate constructor
datas = DataLoader;
% load only folder data
datas = DataLoaderfolder(datas , [1]);
%
session1    =   datas.sessionsDataOffline{1};
% session2    =   datas.getSessionById(2);

 P = session1.P;
%% Data information
 NWindows  = size( P, 1);
 NFreqs    = size( P, 2);
 NChannels = size( P, 3);

 POS = session1.POS;
 DUR = session1.DUR;
 TYP = session1.TYP;
%% Creating vector labels
 CFeedbackPOS =  POS( TYP == 781); % continuos feedback
 CFeedbackDUR =  DUR( TYP == 781); % continuos feedback

 CuePOS =  POS( TYP == 771 |  TYP == 773); % both hands | both feet
 CueDUR =  DUR( TYP == 771 |  TYP == 773); % both hands | both feet
 CueTYP =  TYP( TYP == 771 |  TYP == 773); % both hands | both feet

 FixPOS =  POS( TYP == 786); % Fixation cross
 FixDUR =  DUR( TYP == 786); % Fixation cross
 FixTYP =  TYP( TYP == 786); % fixation cross

 NumTrials = length( CFeedbackPOS);

% We consider the intersting period from Cue apperance to end of continuous feedback
 Ck = zeros( NWindows, 1);
 Tk = zeros( NWindows, 1);
 TrialStart = nan( NumTrials, 1);
 TrialStop  = nan( NumTrials, 1);
 FixStart = nan( NumTrials, 1);
 FixStop  = nan( NumTrials, 1);
for trId = 1: NumTrials
    cstart =  CuePOS(trId);
    cstop  =  CFeedbackPOS(trId) +  CFeedbackDUR(trId) - 1;
     Ck(cstart:cstop) =  CueTYP(trId); % mark TYPE in sample of window 
     Tk(cstart:cstop) = trId; % enumerate sample of window

     TrialStart(trId) = cstart; % mark Pos start
     TrialStop(trId)  = cstop;  % mark Pos stop
     FixStart(trId)   =  FixPOS(trId);
     FixStop(trId)    =  FixPOS(trId) +  FixDUR(trId) - 1;
end

%% Trial extraction

% Extracting data for each trial (be careful that length might be different for few sample)
 MinTrialDur = min( TrialStop -  TrialStart);
 TrialData   = nan( MinTrialDur,  NFreqs,  NChannels,  NumTrials);
 tCk = zeros( NumTrials, 1);
for trId = 1: NumTrials
    cstart =  TrialStart(trId);
    cstop  = cstart +  MinTrialDur - 1;
     TrialData(:, :, :, trId)   =  P(cstart:cstop, :, :);
     tCk(trId) = unique( Ck(cstart:cstop));
end

%% Baseline extraction (from fixation)
 MinFixDur = min( FixStop -  FixStart);
 FixData   = nan( MinFixDur,  NFreqs,  NChannels,  NumTrials);

for trId = 1: NumTrials
    cstart =  FixStart(trId);
    cstop  = cstart +  MinFixDur - 1;
     FixData(:, :, :, trId)   =  P(cstart:cstop, :, :);
end

%% ERD/ERS
% Average and replicate the value of the baseline
 Baseline = repmat(mean( FixData), [size( TrialData, 1) 1 1 1]);
 ERD = log( TrialData./  Baseline);


t = linspace(0, MinTrialDur*datas.wshift, MinTrialDur);

figure;
chandles = [];
ChannelSelected =[7,9,11];
for cId = 1:datas.nclasses
    climits = nan(2, length(ChannelSelected));
    for chId = 1:length(ChannelSelected)
        cdata = mean(ERD(:, :, ChannelSelected(chId), tCk == datas.classId(cId)), 4);
        subplot(2, 3, (cId - 1)*length(ChannelSelected) + chId);
        [chandles,climits(:, chId)] = proc_plotERD_ERS( cdata, datas.channelLb{ChannelSelected(chId)}, t, session1.freqs, datas.classLb{cId});
    end
    
end
set(chandles, 'CLim', [min(min(climits)) max(max(climits))]);

sessionOffline1 = datas.sessionsDataOffline{1};
%% Apply log to the data
SelFreqs = 4:2:48;
fullFreqs = sessionOffline1.freqs;
[freqs, idfreqs] = intersect(fullFreqs, SelFreqs);

U = log(sessionOffline1.P);

NumWins  = size(U, 1);
NumFreqs = size(U, 2);
NumChans = size(U, 3);


%% from ex5_classification_train.m

Rk = sessionOffline1.RkP;
Runs = unique(Rk);
%NumRuns = offlineRuns1;
NumRuns = length(Runs);

%% Labeling the data
disp('[proc] + Labeling the data');

CFeedbackPOS = sessionOffline1.POS(sessionOffline1.TYP == 781);
CFeedbackDUR = sessionOffline1.DUR(sessionOffline1.TYP == 781);

CuePOS = sessionOffline1.POS(sessionOffline1.TYP == 771 | sessionOffline1.TYP == 773 );
CueDUR = sessionOffline1.DUR(sessionOffline1.TYP == 771 | sessionOffline1.TYP == 773 );
CueTYP = sessionOffline1.TYP(sessionOffline1.TYP == 771 | sessionOffline1.TYP == 773 );

NumTrials = length(CueTYP);
% 
% We consider the intersting period from Cue apperance to end of continuous feedback
Ck = zeros(NumWins, 1);
Tk = zeros(NumWins, 1);

for trId = 1:NumTrials
    cstart = CuePOS(trId);
    cstop  = CFeedbackPOS(trId) + CFeedbackDUR(trId) - 1;
    Ck(cstart:cstop) = CueTYP(trId);
    Tk(cstart:cstop) = trId;
end

%% Computing fisher score (for each run)
disp('[proc] + Computing fisher score');
classId = datas.classId(:);
runs = Runs(:);
% disp(classId);
fisherScore = FisherScore(U,Runs,Rk,Ck,classId);
%     score = nan(NumFreqs, NumChans, NumRuns);    
%     NumClasses = length( classId );
%     NumRuns = length(Runs);
%     FS2 = nan(NumFreqs*NumChans, NumRuns);
%     for rId = 1:NumRuns
%         rindex = Rk == Runs(rId); 
%         cmu    = nan(NumFreqs, NumChans, 2);
%         csigma = nan(NumFreqs, NumChans, 2);
%         for cId = 1:NumClasses
%             cindex = rindex & Ck == classId(cId);
%             cmu(:, :, cId) = squeeze(mean(U(cindex, :, :)));
%             csigma(:, :, cId) = squeeze(std(U(cindex, :, :)));
%         end
%         score(:, :, rId) = abs(cmu(:, :, 2) - cmu(:, :, 1)) ./ sqrt( ( csigma(:, :, 1).^2 + csigma(:, :, 2).^2 ) );
%     end
% end

%% Visualization Fisher score
disp('[proc] |- Visualizing fisher score');
OfflineRuns = 1:NumRuns;
climits = [];
handles = nan(NumRuns, 1);
fig1 = figure;
SelChans={};
SelFreqs=[];
FisherScoretemp=fisherScore;
for rId = 1:length(OfflineRuns)
    subplot(1, length(OfflineRuns), rId);

    plotFisherScore(fisherScore, NumRuns,freqs,datas.channelLb);
    % imagesc(fisherScore(:, :, NumRun)');
    % axis square;
    % set(gca, 'XTick', 1:NumFreqs);
    % set(gca, 'XTickLabel', freqs);
    % set(gca, 'YTick', 1:NumChans);
    % set(gca, 'YTickLabel', datas.channelLb);
    % xtickangle(-90);
    
    % title(['Calibration run ' num2str(OfflineRuns(rId))]);
    
    climits = cat(2, climits, get(gca, 'CLim'));
    handles(OfflineRuns(rId)) = gca;
end

% To select the freq and chan with the highest fisher score
for rId = 1:unique(Runs)
    A=FisherScoretemp(:,:,OfflineRuns(rId));
    val = 10;
    while(val>0.9)
        
        [val,idx] = max(A(:));
        if val>0.9
            [row,col] = ind2sub(size(A),idx);
            A(row,col)=0;
            FisherScoretemp(row,col,:)=0;
            SelChans=cat(2,SelChans,datas.channelLb(col));
            SelFreqs=cat(2,SelFreqs,freqs(row));
        end
    end
end


set(handles, 'CLim', [min(min(climits)) max(max(climits))]);

sgtitle('Fisher score');

disp('Topoplot')
load('./Util/chanlocs16.mat');
VisFreq = [11, 12];
figure;
topoplot(mean( fisherScore( 11:12, :,1), 1), chanlocs16, 'headrad', 'rim');
colorbar;
axis image;
title(['Frequency band ' num2str(freqs(VisFreq(1))) '-' num2str(freqs(VisFreq(2))) 'Hz']);

%% Features selection

disp('[proc] |- Select features');
% fprintf('Features: %s \nChannels: %d ',Selchans,SelFreqs);
disp('Features: ');
disp(SelChans);
disp('Channels: ');
disp(SelFreqs);
%SelChans = {'C4', 'C4', 'FC2'};
%SelFreqs = [20 22 22];
% SelChans = unique(SelChans);
NumSelFeatures = length(SelChans);

[~, SelChansId] = ismember(SelChans, datas.channelLb);
[~, SelFreqsId] = ismember(SelFreqs, freqs);
% SelFreqsId=SelFreqs/2;
F = nan(NumWins, NumSelFeatures);
for ftId = 1:NumSelFeatures
    cfrq  = SelFreqsId(ftId);
    cchan = SelChansId(ftId);
    F(:, ftId) = U(:, cfrq, cchan);
end

%% Classifier Train (LDA or QDA);
disp('[proc] + Train classifier');
LabelIdx = Ck == 771 | Ck == 773;
% Model = fitcdiscr(F(LabelIdx, :), Ck(LabelIdx));
Model = fitcdiscr(F(LabelIdx, :), Ck(LabelIdx), 'DiscrimType','quadratic');

%% Classifier accuracy on trainset

[Gk, pp] = predict(Model, F);

SSAcc = 100*sum(Gk(LabelIdx) == Ck(LabelIdx))./length(Gk(LabelIdx));
disp('Model Accuracy:');
disp( SSAcc);

NumClasses = length(classId);
SSClAcc = nan(NumClasses, 1);
for cId = 1:NumClasses
    cindex = Ck == datas.classId(cId);
    SSClAcc(cId) = 100*sum(Gk(cindex) == Ck(cindex))./length(Gk(cindex));
    disp('Model class Accuracy:');
    disp(SSClAcc(cId));
end

%% Saving classifier
disp('[out] + Save classifier');
filename = 'ah7_20201215_classifier.mat';
save(filename, 'Model', 'SelChansId', 'SelFreqsId');

%% Visualize classifier
fig2 = figure;
subplot(1,3,1);
visualizeMu = true;
choose = [1,3];
plot_Classifier(Model,F,LabelIdx,Ck,SelChans,SelFreqs,choose,true);
subplot(1,3,2);
choose = [2,3];
plot_Classifier(Model,F,LabelIdx,Ck,SelChans,SelFreqs,choose, visualizeMu);
subplot(1,3,3);
choose = [1,2];
plot_Classifier(Model,F,LabelIdx,Ck,SelChans,SelFreqs,choose, visualizeMu);


% fig2 = figure;
% h1 = gscatter(F(LabelIdx, 1),F(LabelIdx, 2),Ck(LabelIdx),'rb','ov',[],'off');
% grid on;
% xlim([-8 0]);
% ylim([-8 1.5]);
% xlabel([SelChans{1} '@' num2str(SelFreqs(1)) 'Hz']);
% ylabel([SelChans{2} '@' num2str(SelFreqs(2)) 'Hz']);
% axis square;
% hold on

%% Evidence accumulation
disp('[proc] + Evidence accumulation (exponential smoothing)');
TrialStart = POS(TYP == 781);
NumSamples = size(pp, 1);

ipp = 0.5*ones(size(pp, 1), 1);
alpha = 0.97;

for sId = 2:NumSamples
    
    curr_pp  = pp(sId, 1);
    prev_ipp = ipp(sId-1);
    
    if ismember(sId, TrialStart)
        ipp(sId) = 1./NumClasses;
    else
        ipp(sId) = prev_ipp.*alpha + curr_pp.*(1-alpha);
    end
end


%% Plot accumulated evidence and raw probabilities
fig1 = figure;

% CueClasses    = [771 783 773];
% LineColors = {'b', 'g', 'r'};
% LbClasses     = {'both feet', 'rest', 'both hands'};
% ValueClasses  = [1 0.5 0];
Threshold     = 0.7;

SelTrial = 50;

% Trial 15: good rest
% Trial 80: bad rest
% Trial 55: good both hands
% Trial 58: good both feet
subplot(1,2,1);
plotEvidenceAccumulation(SelTrial, Tk, Ck, pp, ipp, CueTYP, NumTrials, Threshold );
subplot(1,2,2);
plotEvidenceAccumulation(SelTrial, Tk, Ck, pp, ipp, CueTYP, NumTrials, Threshold );
% cindex = Tk == SelTrial;
% [~, ClassIdx] = ismember(unique(Ck(cindex)), CueClasses);
% GreyColor = [150 150 150]/255;
% LineColors = {'b', 'g', 'r'};
% hold on;
% % Plotting raw probabilities
% plot(pp(cindex, 1), 'o', 'Color', GreyColor);
% % Plotting accumulutated evidence
% plot(ipp(cindex), 'k', 'LineWidth', 2);
% % Plotting actual target class
% yline(ValueClasses(ClassIdx), LineColors{ClassIdx}, 'LineWidth', 5);
% % Plotting 0.5 line
% yline(0.5, '--k');
% % Plotting thresholds
% yline(Threshold, 'k', 'Th_{1}');
% yline(1-Threshold, 'k', 'Th_{2}');
% hold off;
% grid on;
% ylim([0 1]);
% xlim([1 sum(cindex)]);
% legend('raw prob', 'integrated prob');
% ylabel('probability/control')
% xlabel('sample');
% title(['Trial ' num2str(SelTrial) '/' num2str(NumTrials) ' - Class ' LbClasses{ClassIdx} ' (' num2str(CueClasses(ClassIdx)) ')']);

%% Compute performances
ActualClass = TYP(TYP == 771 | TYP == 773 | TYP == 783);
Decision = nan(NumTrials, 1);

for trId = 1:NumTrials
    cstart = CFeedbackPOS(trId);
    cstop  = CFeedbackPOS(trId) + CFeedbackDUR(trId) - 1;
    cipp = ipp(cstart:cstop);
    
    endpoint = find( (cipp >= Threshold) | (cipp <= 1 - Threshold), 1, 'first' );
    
    if(isempty(endpoint))
        Decision(trId) = 783;
        continue;
    end
    
    if(cipp(endpoint) >= Threshold)
        Decision(trId) = 771;
    elseif (cipp(endpoint) <= Threshold)
        Decision(trId) = 773;
    end
end

% Removing Rest trials
ActiveTrials = ActualClass ~= 783;
RestTrials = ActualClass == 783;

PerfActive  = 100 * (sum(ActualClass(ActiveTrials) == Decision(ActiveTrials))./sum(ActiveTrials))
PerfResting = 100 * (sum(ActualClass(RestTrials) == Decision(RestTrials))./sum(RestTrials))

RejTrials = Decision == 783;

PerfActive_rej = 100 * (sum(ActualClass(ActiveTrials & ~RejTrials) == Decision(ActiveTrials & ~RejTrials))./sum(ActiveTrials & ~RejTrials))


