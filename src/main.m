%MAIN

% Adding functions class
addpath("./Function/");

% Load Data
data = DataLoader;


sessionsNames = data.getSessionsNames();
sessionsPath = data.getSessionsPaths();
sessionsData = data.getSessionsData();

session = data.getSessionById(1);

