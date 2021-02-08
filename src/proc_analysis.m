%%Proc Script

clearvars;
% Adding functions class
% addpath('./Function/');

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

proc_plotERD_ERS( datas, [7,9,11], t, tCk, session1.freqs, ERD )

proc_plotERD_ERS( datas, [7,9,11], t, tCk, sessionOffline1.freqs, ERD )
