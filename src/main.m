%% MAIN SCRIPT
clear all; close all; clc;
% Adding functions class
addpath('./Function/');
addpath('./Util/');

%% CONFIG PROJECT
datapath       = "../data/";    % Folder with sessions data
f              = 4:2:48;        % SelFreqs 
% Spectrogram params
ml             = 1;             
wl             = 0.5;
ps             = 0.25;                  
ws             = 0.0625;  
wc             = 'backward';
% Classifier training parameters
sc = {'C4','FC2'};
sf = [22 22];

sc1 = {'C1', 'Cz', 'Cz'};
sf1 = [12 24 22];

%% Create data presenter instance
presenter = DataPresenter();

%% Load Data
data        = DataLoader(datapath,f,ml,wl,ps,ws,wc,presenter);

%% Process Data
processor   = DataProcessing(data,f,presenter);

%% Classifier for data
%classifier  = DataClassifier(processor,sc1,sf1,presenter);
classifier  = DataClassifier(processor,sc,sf,presenter);%,SelChans,channelLb,SelFreqs,freqs,U,Ck,NumWins);



%{
%% EEG processing - Train classifier
%
%% Implementation steps:
%
% 1 + Import PSD data
%
% 2 + Concatenate PSD data
%
% 3 + Perform feature selection
%
% 4 + Train classifier
%
% 5 + Save classifier

%% General information
%% General information
clearvars; clc;
fprintf("Loading Biosig v3.7.2...\n");
addpath('./Util/biosig4octmat-3.7.2/biosig/t200_FileAccess','./Util/biosig4octmat-3.7.2/biosig/t210_Events','./Util/biosig4octmat-3.7.2/biosig/t250_ArtifactPreProcessingQualityControl','./Util/biosig4octmat-3.7.2/biosig/t300_FeatureExtraction','./Util/biosig4octmat-3.7.2/biosig/t310_ERDSMaps','./Util/biosig4octmat-3.7.2/biosig/t320_Nirs','./Util/biosig4octmat-3.7.2/biosig/t330_StimFit','./Util/biosig4octmat-3.7.2/biosig/t400_Classification','./Util/biosig4octmat-3.7.2/biosig/t450_MultipleTestStatistic','./Util/biosig4octmat-3.7.2/biosig/t490_EvaluationCriteria','./Util/biosig4octmat-3.7.2/biosig/t500_Visualization','./Util/biosig4octmat-3.7.2/biosig/t501_VisualizeCoupling');

datapath   = '../data/';
channelLb  = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};
channelId  = 1:length(channelLb);
classId    = [773 771];   
classLb    = {'both hands', 'both feet'};
nclasses   = length(classId);
modalityId = [0 1];
modalityLb = {'offline', 'online'};

mlength    = 1;
wlength    = 0.5;
pshift     = 0.25;                  
wshift     = 0.0625;  
selfreqs   = 4:2:96;
winconv = 'backward'; 


%GDFfile = fullfile(datapath, 'ah7.20170613.161402.offline.mi.mi_bhbf.gdf');
%GDFfile = fullfile(datapath, 'ah7.20170613.162331.offline.mi.mi_bhbf.gdf');
GDFfile = fullfile(datapath, 'ah7.20170613.162934.offline.mi.mi_bhbf.gdf');


%% Loading and concatenate GDF file
disp(['[io] - Loading GDF file : ' GDFfile]);
[s, h] = sload(GDFfile);
s = s(:, channelId);
SampleRate = h.SampleRate;

%% Spatial filters
disp('[proc] |- Applying CAR and Laplacian');
load('laplacian16.mat');
s_lap = s*lap;

%% Spectrogram
disp('[proc] |- Computing spectrogram');
[P, freqgrid] = proc_spectrogram(s_lap, wlength, wshift, pshift, SampleRate, mlength);  

%% Selecting desired frequencies
[freqs, idfreqs] = intersect(freqgrid, selfreqs);
P = P(:, idfreqs, :);
event = h.EVENT;

%% Extracting events
disp('[proc] |- Extract and convert the events');
events.TYP = h.EVENT.TYP;
events.POS = proc_pos2win(h.EVENT.POS, wshift*h.SampleRate, winconv, mlength*h.SampleRate);
events.DUR = floor(h.EVENT.DUR/(wshift*h.SampleRate)) + 1;
events.conversion = winconv;
    
%% Data information
NWindows  = size(P, 1);
NFreqs    = size(P, 2);
NChannels = size(P, 3);

%% Creating vector labels
CFeedbackPOS = events.POS(events.TYP == 781);
CFeedbackDUR = events.DUR(events.TYP == 781);

CuePOS = events.POS(events.TYP == 771 | events.TYP == 773);
CueDUR = events.DUR(events.TYP == 771 | events.TYP == 773);
CueTYP = events.TYP(events.TYP == 771 | events.TYP == 773);

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
for cId = 1:nclasses
    
    climits = nan(2, length(ChannelSelected));
    for chId = 1:length(ChannelSelected)
        subplot(2, 3, (cId - 1)*length(ChannelSelected) + chId);
        cdata = mean(ERD(:, :, ChannelSelected(chId), tCk == classId(cId)), 4);
        imagesc(t, freqs, cdata');
        set(gca,'YDir','normal');
        climits(:, chId) = get(gca, 'CLim');
        chandles = cat(1, chandles, gca);
        colormap(hot);
        colorbar;
        title(['Channel ' channelLb{ChannelSelected(chId)} ' | ' classLb{cId}]);
        xlabel('Time [s]');
        ylabel('Frequency [Hz]');
        line([1 1],get(gca,'YLim'),'Color',[0 0 0])
    end
    
end
set(chandles, 'CLim', [min(min(climits)) max(max(climits))]);

%save('ah7.20170613.161402.offline.mi.mi_bhbf.mat');
%save('ah7.20170613.162331.offline.mi.mi_bhbf.mat');
save('ah7.20170613.162934.offline.mi.mi_bhbf.mat');
%}











%{


%clearvars; clc;

channelLb  = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};
channelId  = 1:length(channelLb);
classId    = [773 771];   
classLb    = {'both hands', 'both feet'};
nclasses   = length(classId);
modalityId = [0 1];
modalityLb = {'offline', 'online'};

%% Importing EEG data
files{1} = 'ah7.20170613.161402.offline.mi.mi_bhbf.mat';
files{2} = 'ah7.20170613.162331.offline.mi.mi_bhbf.mat';
files{3} = 'ah7.20170613.162934.offline.mi.mi_bhbf.mat';

nfiles = length(files);

% Loading and concatenate PSD files and events
P = [];
TYP = [];
POS = [];
DUR = [];
Rk = []; Mk = [];

for fId = 1:nfiles
    disp(['[io] - Loading PSD file ' num2str(fId) '/' num2str(nfiles) ': ' files{fId}]);
    cfilename = files{fId}; 
    cdata = load(cfilename);
    
    % Extract the current PSD
    cpsd = cdata.P;
    
    % Extract the current events
    cevents = cdata.events;
    TYP = cat(1, TYP, cevents.TYP);
    DUR = cat(1, DUR, cevents.DUR);
    POS = cat(1, POS, cevents.POS + size(P, 1));
    
    % Create Rk vector (run)
    cRk = fId*ones(size(cpsd, 1), 1);
    Rk = cat(1, Rk, cRk);
    
    % Create Mk vector (modality)
    if( contains(cfilename, 'offline') == true)
        cMk = modalityId(1)*ones(size(cpsd, 1), 1);
    elseif ( contains(cfilename, 'online') == true )
        cMk = modalityId(2)*ones(size(cpsd, 1), 1);
    else
        error(['Unknown modality for run: ' cfilename]);
    end
    Mk = cat(1, Mk, cMk);
    
    P = cat(1, P, cpsd);
    fullfreqs = cdata.freqs; 
    SampleRate = cdata.SampleRate;
end

%% Apply log to the data

SelFreqs = 4:2:48;
[freqs, idfreqs] = intersect(fullfreqs, SelFreqs);

U = log(P(:, idfreqs, :));

NumWins  = size(U, 1);
NumFreqs = size(U, 2);
NumChans = size(U, 3);

Runs = unique(Rk);
NumRuns = length(Runs);

%% Labeling the data
disp('[proc] + Labeling the data');

CFeedbackPOS = POS(TYP == 781);
CFeedbackDUR = DUR(TYP == 781);

CuePOS = POS(TYP == 771 | TYP == 773 | TYP == 783);
CueDUR = DUR(TYP == 771 | TYP == 773 | TYP == 783);
CueTYP = TYP(TYP == 771 | TYP == 773 | TYP == 783);

NumTrials = length(CueTYP);

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
Classes = [771 773];
NumClasses = length(Classes);

FisherScore = nan(NumFreqs, NumChans, NumRuns);
FS2 = nan(NumFreqs*NumChans, NumRuns);
for rId = 1:NumRuns
    rindex = Rk == Runs(rId); 
    
    cmu    = nan(NumFreqs, NumChans, 2);
    csigma = nan(NumFreqs, NumChans, 2);
    
    for cId = 1:NumClasses
        cindex = rindex & Ck == Classes(cId);
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
    set(gca, 'YTickLabel', channelLb);
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

[~, SelChansId] = ismember(SelChans, channelLb);
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

%[Gk, pp] = predict(Model, F);

%SSAcc = 100*sum(Gk(LabelIdx) == Ck(LabelIdx))./length(Gk(LabelIdx));

%SSClAcc = nan(NumClasses, 1);
%for cId = 1:NumClasses
%    cindex = Ck == Classes(cId);
%    SSClAcc(cId) = 100*sum(Gk(cindex) == Ck(cindex))./length(Gk(cindex));
%end

%% Saving classifier
%disp('[out] + Save classifier');
%filename = 'ah7_20201215_classifier.mat';
%save(filename, 'Model', 'SelChansId', 'SelFreqsId');

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
%}

































%classifier.loadFromProcessor()



%{
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


%sessionOffline1 = data.getSessionOfflineById(9);

%session2        = data.getSessionByName("20190711_F1");         % Load data of session 20190711_F1
%sessionOnline2  = data.getSessionOnlineByName("20190711_F1");   % Load data of session 20190711_F1
%sessionOffline2 = data.getSessionOfflineByName("20190711_F1");  % Load data of session 20190711_F1



%% Example of data usage from a sessions

disp(session1.SampleRate)     % Display Sample rate
figure;
subplot(3,2,1);
s1 = session1.s(:,1:16); % because the column 17 is empty
plot(sessionOffline1.P);   % Plot samples
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




%% Spatial filters
disp('[proc] |- Applying CAR and Laplacian');
load('./Util/laplacian16.mat');
%s1 = session1.s(:,1:16); % because the column 17 is empty
s_lap = s1*lap;


%% Spectrogram (PSD)
disp('[proc] |- Computing spectrogram');
[P, freqgrid] = proc_spectrogram(s_lap, wlength, wshift, pshift, session1.SampleRate, mlength);  
size(P)

%% Selecting desired frequencies
[freqs, idfreqs] = intersect(freqgrid, selfreqs);
P = P(:, idfreqs, :);
size(P)
%% Extracting events
disp('[proc] |- Extract and convert the events');
events.TYP = session1.TYP;
events.POS = proc_pos2win(session1.POS, wshift*session1.SampleRate, winconv, mlength*session1.SampleRate);
events.DUR = floor(session1.DUR/(wshift*session1.SampleRate)) + 1;
events.conversion = winconv;

%% Data information
P = sessionOffline1.P;
NWindows  = size(P, 1);
NFreqs    = size(P, 2);
NChannels = size(P, 3);

%% Creating vector labels
CFeedbackPOS = sessionOffline1.POS(sessionOffline1.TYP == 781);
CFeedbackDUR = sessionOffline1.DUR(sessionOffline1.TYP == 781);

CuePOS = sessionOffline1.POS(sessionOffline1.TYP == 771 | sessionOffline1.TYP == 773 );
CueDUR = sessionOffline1.DUR(sessionOffline1.TYP == 771 | sessionOffline1.TYP == 773 );
CueTYP = sessionOffline1.TYP(sessionOffline1.TYP == 771 | sessionOffline1.TYP == 773);

FixPOS = sessionOffline1.POS(sessionOffline1.TYP == 786);
FixDUR = sessionOffline1.DUR(sessionOffline1.TYP == 786);
FixTYP = sessionOffline1.TYP(sessionOffline1.TYP == 786);

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
        imagesc(t, session1.freqs, cdata');
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
OfflineRuns = 1:NumRuns;
climits = [];
handles = nan(NumRuns, 1);
fig1 = figure;
SelChans={};
SelFreqs=[];
FisherScoretemp=FisherScore;
for rId = 1:length(OfflineRuns)
    subplot(1, length(OfflineRuns), rId);
    imagesc(FisherScore(:, :, OfflineRuns(rId))');
    
    % To select the freq and chan with the highest fisher score
    A=FisherScoretemp(:,:,OfflineRuns(rId));
    val = 10;
    while(val>0.9)
        
        [val,idx] = max(A(:));
        if val>0.9
            [row,col] = ind2sub(size(A),idx);
            A(row,col)=0;
            FisherScoretemp(row,col,:)=0;
            SelChans=cat(2,SelChans,data.channelLb(col));
            SelFreqs=cat(2,SelFreqs,freqs(row));
        end
    end
    
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
%SelChans = {'C4', 'C4', 'FC2'};
%SelFreqs = [20 22 22];

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





%%%%%%
%% from ex6_control_framework_intergartion.m

%% Apply log to the data

SelFreqs = 4:2:48;
sessionOnline = data.getSessionOnlineById(1);   
fullfreqs = sessionOnline.freqs;
[freqs, idfreqs] = intersect(fullfreqs, SelFreqs);

U = log(sessionOnline.P(:, idfreqs, :));

NumWins  = size(U, 1);
NumFreqs = size(U, 2);
NumChans = size(U, 3);

Runs = unique(sessionOnline.RkP);
NumRuns = length(Runs);

TYP= sessionOnline.TYP;
POS=sessionOnline.POS;
DUR=sessionOnline.DUR;

%% Labeling the data
disp('[proc] + Labeling the data');

CFeedbackPOS = POS(TYP == 781);
CFeedbackDUR = DUR(TYP == 781);

CuePOS = POS(TYP == 771 | TYP == 773 | TYP == 783);
CueDUR = DUR(TYP == 771 | TYP == 773 | TYP == 783);
CueTYP = TYP(TYP == 771 | TYP == 773 | TYP == 783);

NumTrials = length(CueTYP);

% We consider the intersting period from Cue apperance to end of continuous feedback
Ck = zeros(NumWins, 1);
Tk = zeros(NumWins, 1);
for trId = 1:NumTrials
    cstart = CuePOS(trId);
    cstop  = CFeedbackPOS(trId) + CFeedbackDUR(trId) - 1;
    Ck(cstart:cstop) = CueTYP(trId);
    
    cstart_cfb = CFeedbackPOS(trId);
    Tk(cstart_cfb:cstop) = trId;
end

%% Loading the classifier
classifier_name = 'ah7_20201215_classifier.mat';
MDL = load(classifier_name);
SelChansId = MDL.SelChansId;
SelFreqsId = MDL.SelFreqsId;
Model = MDL.Model;

%% Features extraction
disp('[proc] |- Extracing features');

NumSelFeatures = length(SelChansId);

F = nan(NumWins, NumSelFeatures);
for ftId = 1:NumSelFeatures
    cfrq  = SelFreqsId(ftId);
    cchan = SelChansId(ftId);
    F(:, ftId) = U(:, cfrq, cchan);
end

%% Classifier Test
disp('[proc] + Evaluate classifier');
LabelIdx = Ck == 771 | Ck == 773;

%% Overall classifier accuracy on testset

[Gk, pp] = predict(Model, F);

SSAcc = 100*sum(Gk(LabelIdx) == Ck(LabelIdx))./length(Gk(LabelIdx));
Classes = classId;
NumClasses = length(classId);
SSClAcc = nan(NumClasses, 1);
for cId = 1:NumClasses
    cindex = Ck == Classes(cId);
    SSClAcc(cId) = 100*sum(Gk(cindex) == Ck(cindex))./length(Gk(cindex));
end

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

CueClasses    = [771 783 773];
LbClasses     = {'both feet', 'rest', 'both hands'};
ValueClasses  = [1 0.5 0];
Threshold     = 0.7;

SelTrial = 50;

% Trial 15: good rest
% Trial 80: bad rest
% Trial 55: good both hands
% Trial 58: good both feet

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
title(['Trial ' num2str(SelTrial) '/' num2str(NumTrials) ' - Class ' LbClasses{ClassIdx} ' (' num2str(CueClasses(ClassIdx)) ')']);


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

%}
