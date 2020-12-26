%% MAIN SCRIPT

% Adding functions class
addpath("./Function/");

% Load Data
data = DataLoader;


%% Example how to get sessions data and names from loader

sessionsNames = data.getSessionsNames();  % Load Sessions names
sessionsPath  = data.getSessionsPaths();  % Load Sessions Path
sessionsData  = data.getSessionsData();   % Load Sessions Data

session1      = data.getSessionById(1);                 % Load data of session #1
session2      = data.getSessionByName("20190711_F1");   % Load data of session 20190711_F1

%% Example of data usage from a sessions
%{
disp(session1.SampleRate)     % Display Sample rate

plot(session1.s);             % Plot samples
plot(session1.TYP);           % Plot session TYP vector 
plot(session1.DUR);           % Plot session DUR vector
plot(session1.POS);           % Plot session POS vector 
plot(session1.Rk);            % Plot session Rk 
plot(session1.Mk);            % Plot session Mk 
%}
