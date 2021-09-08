%% MAIN SCRIPT
clear all; close all; clc;  % Clear prev environment
addpath('./Function/');     % Add functions path
addpath('./Util/');         % Add util path with biosig and laplacian

%% CONFIG PROJECT

%% Configure input
datapath       = "../data/";    % Folder with sessions data (each session must be in subfolder that contains related gdf files)
f              = 4:2:48;        % SelFreqs [ 4 ~ 48 with step 2]

%% Spectrogram params
ml             = 1;             
wl             = 0.5;
ps             = 0.25;                  
ws             = 0.0625;  
wc             = 'backward';

%% Classifier training parameters
sc = {'C4','FC2'};
sf = [22 22];

%% output configuration
display_input_data           = 0; % (0 = Not Display | 1 = Display if possible)
display_erd_ers              = 0; % (0 = Not Display | 1 = Display if possible)
display_fisher_score         = 0; % (0 = Not Display | 1 = Display if possible)
display_classifier           = 0; % (0 = Not Display | 1 = Display if possible)
display_accumulated_evidence = 1; % (0 = Not Display | 1 = Display if possible)

%% PROJECT RUN


%% Create data presenter instance
presenter   = DataPresenter(display_input_data,display_erd_ers,display_fisher_score,display_classifier,display_accumulated_evidence);

%% Load Data
data        = DataLoader(datapath,f,ml,wl,ps,ws,wc,presenter);

%% Process Data
processor   = DataProcessing(data,f,presenter);

%% Classifier for data
classifier  = DataClassifier(processor,sc,sf,presenter);

