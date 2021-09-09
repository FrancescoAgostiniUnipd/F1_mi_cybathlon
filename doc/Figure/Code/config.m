%% CONFIG PROJECT
datapath       = "../data/";    % Folder with sessions data
f              = 4:2:48;        % SelFreqs 
% Spectrogram params
ml             = 1;             % mlength        
wl             = 0.5;           % wlength
ps             = 0.25;          % pshift         
ws             = 0.0625;        % wshift
wc             = 'backward';    % winconv
% Feature selection params
sc             = {'C4','FC2'};
sf             = [22 22];
% output configuration > (0 = Not Display | 1 = Display if possible)
display_input_data           = 0;  
display_erd_ers              = 0; 
display_fisher_score         = 0; 
display_classifier           = 0; 
display_accumulated_evidence = 0; 