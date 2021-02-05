%% MAIN SCRIPT

% Adding functions class
addpath('./Function/');
addpath('./Util/');

% Load Data
data = DataLoader;


%% Example how to get sessions data and names from loader

sessionsNames = data.getSessionsNames();  % Load Sessions names
sessionsPath  = data.getSessionsPaths();  % Load Sessions Path
sessionsData  = data.getSessionsData();   % Load Sessions Data

% Session 1 example
session1        = data.getSessionById(1);                 % Load all data of session #1
sessionOnline1  = data.getSessionOnlineById(1);           % Load all online data of session #1
sessionOffline1 = data.getSessionOfflineById(1);          % Load all offline data of session #1
allRuns1        = data.allRuns{1};                        % overall runs in session 1
offlineRuns1    = data.offlineRuns{1};                    % offline runs in session 1
onlineRuns1     = data.onlineRuns{1};                     % online runs in session 1


session2        = data.getSessionByName("20190711_F1");         % Load data of session 20190711_F1
sessionOnline2  = data.getSessionOnlineByName("20190711_F1");   % Load data of session 20190711_F1
sessionOffline2 = data.getSessionOfflineByName("20190711_F1");  % Load data of session 20190711_F1



%% Example of data usage from a sessions

disp(session1.SampleRate)     % Display Sample rate
figure;
subplot(3,2,1);
s1 = session1.s(:,1:16); % because the column 17 is empty
plot(s1);   % Plot samples
title('samples')
subplot(3,2,2);
plot(session1.TYP);           % Plot session TYP vector 
title('session TYP')
subplot(3,2,3);
plot(session1.DUR);           % Plot session DUR vector
title('session DUR')
subplot(3,2,4);
plot(session1.POS);           % Plot session POS vector 
title('session POS')
subplot(3,2,5);
plot(session1.Rk);            % Plot session Rk 
title('session Rk')
subplot(3,2,6);
plot(session1.Mk);            % Plot session Mk 
title('session Mk')


%% ????????? Data analize ??????  from ex3_spectogram_students.com
mlength    = 1; % ??
wlength    = 0.5; % ??
pshift     = 0.25;  % ??              
wshift     = 0.0625; % ??
winconv = 'backward'; % ??
selfreqs   = 4:2:96; % ??

%% Spatial filters
disp('[proc] |- Applying CAR and Laplacian');
load('./Util/laplacian16.mat');
%s1 = session1.s(:,1:16); % because the column 17 is empty
s_lap = s1*lap;



%% Spectrogram (PSD)
disp('[proc] |- Computing spectrogram');
[P, freqgrid] = proc_spectrogram(s_lap, wlength, wshift, pshift, session1.SampleRate, mlength);  


%% Selecting desired frequencies
[freqs, idfreqs] = intersect(freqgrid, selfreqs);
P = P(:, idfreqs, :);

%% Extracting events
disp('[proc] |- Extract and convert the events');
events.TYP = session1.TYP;
events.POS = proc_pos2win(session1.POS, wshift*session1.SampleRate, winconv, mlength*session1.SampleRate);
events.DUR = floor(session1.DUR/(wshift*session1.SampleRate)) + 1;
events.conversion = winconv;

%% Data information
NWindows  = size(P, 1);
NFreqs    = size(P, 2);
NChannels = size(P, 3);

%% Creating vector labels
CFeedbackPOS = events.POS(events.TYP == 781);
CFeedbackDUR = events.DUR(events.TYP == 781);

CuePOS = events.POS(events.TYP == 771 | events.TYP == 773 | events.TYP == 783);
CueDUR = events.DUR(events.TYP == 771 | events.TYP == 773 | events.TYP == 783);
CueTYP = events.TYP(events.TYP == 771 | events.TYP == 773 | events.TYP == 783);

FixPOS = events.POS(events.TYP == 786);
FixDUR = events.DUR(events.TYP == 786);
FixTYP = events.TYP(events.TYP == 786);

NumTrials = length(CFeedbackPOS);

% We consider the intersting period from Cue apperance to end of continuous feedback
Ck = zeros(NWindows, 1);
Tk = zeros(NWindows, 1);
TrialStart = nan(NumTrials, 1);
TrialStop  = nan(NumTrials, 1);
FixStart = nan(NumTrials, 1);
FixStop  = nan(NumTrials, 1);
for trId = 1:NumTrials
    cstart = CuePOS(trId);
    cstop  = CFeedbackPOS(trId) + CFeedbackDUR(trId) - 1;
    Ck(cstart:cstop) = CueTYP(trId);
    Tk(cstart:cstop) = trId;
    
    TrialStart(trId) = cstart;
    TrialStop(trId)  = cstop;
    FixStart(trId)   = FixPOS(trId);
    FixStop(trId)    = FixPOS(trId) + FixDUR(trId) - 1;
end

%% Trial extraction

% Extracting data for each trial (be careful that length might be different for few sample)
disp('[proc] + Extracting data for each trial');
MinTrialDur = min(TrialStop - TrialStart);
TrialData   = nan(MinTrialDur, NFreqs, NChannels, NumTrials);
tCk = zeros(NumTrials, 1);
for trId = 1:NumTrials
    cstart = TrialStart(trId);
    cstop  = cstart + MinTrialDur - 1;
    TrialData(:, :, :, trId)   = P(cstart:cstop, :, :);
   
    tCk(trId) = unique(Ck(cstart:cstop));
end

%% Baseline extraction (from fixation)
disp('[proc] + Extracting baseline data for each trial');
MinFixDur = min(FixStop - FixStart);
FixData   = nan(MinFixDur, NFreqs, NChannels, NumTrials);

for trId = 1:NumTrials
    cstart = FixStart(trId);
    cstop  = cstart + MinFixDur - 1;
    FixData(:, :, :, trId)   = P(cstart:cstop, :, :);
end

%% ERD/ERS
disp('[proc] + Computing ERD/ERS');

% Average and replicate the value of the baseline
Baseline = repmat(mean(FixData), [size(TrialData, 1) 1 1 1]);
ERD = log(TrialData./ Baseline);

%% Visualization 1
figure;
t = linspace(0, MinTrialDur*wshift, MinTrialDur);
ChannelSelected = [7 9 11]; 

chandles = [];
for cId = 1:data.nclasses
    
    climits = nan(2, length(ChannelSelected));
    for chId = 1:length(ChannelSelected)
        subplot(2, 3, (cId - 1)*length(ChannelSelected) + chId);
        cdata = mean(ERD(:, :, ChannelSelected(chId), tCk == data.classId(cId)), 4);
        imagesc(t, freqs, cdata');
        set(gca,'YDir','normal');
        climits(:, chId) = get(gca, 'CLim');
        chandles = cat(1, chandles, gca);
        colormap(hot);
        colorbar;
        title(['Channel ' data.channelLb{ChannelSelected(chId)} ' | ' data.classLb{cId}]);
        xlabel('Time [s]');
        ylabel('Frequency [Hz]');
        line([1 1],get(gca,'YLim'),'Color',[0 0 0])
    end
    
end
set(chandles, 'CLim', [min(min(climits)) max(max(climits))]);



%% from ex4_features_selection_students.m

%% Apply log to the data
SelFreqs = 4:2:48;
fullFreqs = sessionOffline1.freqs;
[freqs, idfreqs] = intersect(fullfreqs, SelFreqs);

U = log(P);

NumWins  = size(U, 1);
NumFreqs = size(U, 2);
NumChans = size(U, 3);


%% from ex5_classification_train.m
%Rk = ones(size(U, 1), 1);
Rk = sessionOffline1.Rk;
Runs = unique(Rk);
NumRuns = offlineRuns1;
%NumRuns = length(Runs);

% %% Labeling the data
% disp('[proc] + Labeling the data');
% 
% CFeedbackPOS = events.POS(events.TYP == 781);
% CFeedbackDUR = events.DUR(events.TYP == 781);
% 
% CuePOS = events.POS(events.TYP == 771 | events.TYP == 773 | events.TYP == 783);
% CueDUR = events.DUR(events.TYP == 771 | events.TYP == 773 | events.TYP == 783);
% CueTYP = events.TYP(events.TYP == 771 | events.TYP == 773 | events.TYP == 783);
% 
% NumTrials = length(CueTYP);
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

NumClasses = length(data.classId);

FisherScore = nan(NumFreqs, NumChans, NumRuns);
FS2 = nan(NumFreqs*NumChans, NumRuns);
for rId = 1:NumRuns
    rindex = Rk == Runs(rId); 
    
    cmu    = nan(NumFreqs, NumChans, 2);
    csigma = nan(NumFreqs, NumChans, 2);
    
    for cId = 1:NumClasses
        cindex = rindex & Ck == data.classId(cId);
        cmu(:, :, cId) = squeeze(mean(U(cindex, :, :)));
        csigma(:, :, cId) = squeeze(std(U(cindex, :, :)));
    end
    
    FisherScore(:, :, rId) = abs(cmu(:, :, 2) - cmu(:, :, 1)) ./ sqrt( ( csigma(:, :, 1).^2 + csigma(:, :, 2).^2 ) );
end

%% Visualization Fisher score
disp('[proc] |- Visualizing fisher score');
OfflineRuns = [1 2 3];
climits = [];
handles = nan(NumRuns, 1);
fig1 = figure;

for rId = 1:length(OfflineRuns)
    subplot(1, 3, rId);
    imagesc(FisherScore(:, :, OfflineRuns(rId))');
    axis square;
    set(gca, 'XTick', 1:NumFreqs);
    set(gca, 'XTickLabel', freqs);
    set(gca, 'YTick', 1:NumChans);
    set(gca, 'YTickLabel', data.channelLb);
    xtickangle(-90);
    
    title(['Calibration run ' num2str(OfflineRuns(rId))]);
    
    climits = cat(2, climits, get(gca, 'CLim'));
    handles(OfflineRuns(rId)) = gca;
end


set(handles, 'CLim', [min(min(climits)) max(max(climits))]);

sgtitle('Fisher score');

%% Features selection
disp('[proc] |- Select features');
SelChans = {'C1', 'Cz', 'Cz'};
SelFreqs = [12 24 22];
NumSelFeatures = length(SelChans);

[~, SelChansId] = ismember(SelChans, data.channelLb);
[~, SelFreqsId] = ismember(SelFreqs, freqs);

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

SSClAcc = nan(NumClasses, 1);
for cId = 1:NumClasses
    cindex = Ck == data.classId(cId);
    SSClAcc(cId) = 100*sum(Gk(cindex) == Ck(cindex))./length(Gk(cindex));
end

%% Saving classifier
disp('[out] + Save classifier');
filename = 'ah7_20201215_classifier.mat';
save(filename, 'Model', 'SelChansId', 'SelFreqsId');

%% Visualize classifier
fig2 = figure;
h1 = gscatter(F(LabelIdx, 1),F(LabelIdx, 2),Ck(LabelIdx),'kb','ov^',[],'off');
grid on;
xlim([-8 0]);
ylim([-8 1.5]);
xlabel([SelChans{1} '@' num2str(SelFreqs(1)) 'Hz']);
ylabel([SelChans{2} '@' num2str(SelFreqs(2)) 'Hz']);
axis square;
hold on

% Linear
% K = Model.Coeffs(1,2).Const;
% L = Model.Coeffs(1,2).Linear;
% f = @(x1,x2) K + L(1)*x1 + L(2)*x2;

% Quadratic
K = Model.Coeffs(1,2).Const;
L = Model.Coeffs(1,2).Linear;
Q = Model.Coeffs(1,2).Quadratic;
f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;

h2 = fimplicit(f);
h2.Color = 'r';
h2.LineWidth = 2;
h2.DisplayName = 'Boundary between boht hands & both feet';
legend('both feet', 'both hands', 'Boundary');
hold off;

