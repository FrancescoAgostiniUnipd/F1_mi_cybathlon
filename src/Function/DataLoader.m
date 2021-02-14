%% DATA LOADER FUNCTION CLASS 
%{
Description: Collection of utils and functions userd for load 
              and chain .gdf EEG description file on workspace

Authors: Agostini Francesco (francesco.agostini.5@studenti.unipd.it)
          Deschaux Oph�lie   (opheliecandicemarine.deschaux@studenti.unipd.it)
          Marcon Francesco   (francesco.marcon.2@studenti.unipd.it)

Version: 0.4

%}
%% Class Dataloader definition
classdef DataLoader
    %% Class Properties 
    properties
        datapath               % Path where data are stored
        datasample             % Sample frequency
        channelLb              % Label fot EEG channels
        channelId              % Channel ID number
        classId                % ID of classes vector
        classLb                % Classes description vector
        classAllId             % ID of classes vector ALL
        classAllLb             % Classes description vector ALL
        nclasses               % N of classes
        modalityId             % ID of modality
        modalityLb             % Label of modality (online ~ offline)
        mlength
        wlength 
        pshift                    
        wshift    
        selfreqs
        winconv 
        sessionsNames          % Sessions names present in data directory
        sessionsPaths          % sessions path correspondent to sessions names
        sessionsData           % sessions data correspondent to sessions names
        sessionsDataOffline    % offline sessions data correspondent to sessions names
        sessionsDataOnline     % online sessions data correspondent to sessions names
        allRuns                % overall runs in session
        offlineRuns            % offline runs in session
        onlineRuns             % online runs in session
    end
    
    %% Class constructor and loader method
    methods
        %% Constructor
        function obj = DataLoader()
            % Setting up env
            % obj.datapath            = '../datalectures/';
            obj.datapath            = '../data/';
            obj.datasample          = 512;
            obj.channelLb           = {'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'};
            obj.channelId           = 1:length(obj.channelLb);
            obj.classAllId          = [786 773 771 783 781 897 898];
            obj.classId             = [773 771];
            obj.classAllLb          = {'Fixation cross','Both hands','Both feet','Rest','Continuos feedback','Target hint','Target miss'}; % we have to add class labels
            obj.classLb             = {'Both hands','Both feet'}; % we have to add class labels
            obj.nclasses            = length(obj.classId);
            obj.modalityId          = [0 1];
            obj.modalityLb          = {'offline','online'};
            obj.mlength             = 1;
            obj.wlength             = 0.5;
            obj.pshift              = 0.25;                  
            obj.wshift              = 0.0625;  
            obj.selfreqs            = 4:2:48;
            obj.winconv             = 'backward'; 
            obj.sessionsData        = [];
            obj.sessionsDataOffline = [];
            obj.sessionsDataOnline  = [];
            obj.allRuns             = [];
            obj.offlineRuns         = [];
            obj.onlineRuns          = [];
            
            % Start loading data names
            sprintf("\n\nStart loading sessions ... \n");
            obj = obj.listSessions();
            % Start loading gdf files
            % obj = obj.loadSessions();
            
            
            sprintf("\nAll %d sessions loaded correctly!\n\n",length(obj.sessionsNames));
        end

        %% Load only folder number {selectedFolderArray}

        % for use comment loadSessions() in the constructor DataLoader()
        
        %     obj: is costructor object
        %     folder: is a array of select folder
        function obj = DataLoaderfolder( obj, folder)
            obj = DataLoader();
            
            obj.loadBiosig();

            for i = folder
                fprintf("Loading session #%d %s \n",i,obj.sessionsPaths{i});
                obj = obj.loadSession(obj.sessionsPaths{i},i);
            end
            fprintf("\nAll %d sessions loaded correctly!\n\n",length(obj.sessionsData));
        end
        
        
        function [curr_P,curr_freqs,curr_TYP,curr_DUR,curr_POS] = preprocessing(obj,curr_s,curr_h)
            curr_SampleRate = curr_h.SampleRate;
            s = curr_s(:, 1:16);

            % Applying CAR and Laplacian
            load('./Util/laplacian16.mat');
            curr_s_lap = s*lap;
            
            % Computing spectrogram            
            [curr_P, curr_freqgrid] = proc_spectrogram(curr_s_lap, obj.wlength, obj.wshift, obj.pshift, curr_SampleRate, obj.mlength);
            
            %% Selecting desired frequencies
            [curr_freqs, curr_idfreqs] = intersect(curr_freqgrid, obj.selfreqs);
            curr_P = curr_P(:, curr_idfreqs, :);
            
            %% Extracting events
            curr_TYP = curr_h.EVENT.TYP;
            curr_POS = proc_pos2win(curr_h.EVENT.POS, obj.wshift*curr_h.SampleRate, obj.winconv, obj.mlength*curr_h.SampleRate);
            curr_DUR = floor(curr_h.EVENT.DUR/(obj.wshift*curr_h.SampleRate)) + 1;
            events.conversion = obj.winconv;
        end
        
        function [curr_P,curr_freqs,curr_ERD,curr_TYP,curr_DUR,curr_POS] = ERDS(obj,curr_s,curr_h)
%             curr_SampleRate = curr_h.SampleRate;
%             s = curr_s(:, 1:16);
% 
%             % Applying CAR and Laplacian
%             load('./Util/laplacian16.mat');
%             curr_s_lap = s*lap;
%             
%             % Computing spectrogram            
%             [curr_P, curr_freqgrid] = proc_spectrogram(curr_s_lap, obj.wlength, obj.wshift, obj.pshift, curr_SampleRate, obj.mlength);
%             
%             %% Selecting desired frequencies
%             [curr_freqs, curr_idfreqs] = intersect(curr_freqgrid, obj.selfreqs);
%             curr_P = curr_P(:, curr_idfreqs, :);
%             
%             %% Extracting events
%             curr_TYP = curr_h.EVENT.TYP;
%             curr_POS = proc_pos2win(curr_h.EVENT.POS, obj.wshift*curr_h.SampleRate, obj.winconv, obj.mlength*curr_h.SampleRate);
%             curr_DUR = floor(curr_h.EVENT.DUR/(obj.wshift*curr_h.SampleRate)) + 1;
%             events.conversion = obj.winconv;


            [curr_P,curr_freqs,curr_TYP,curr_DUR,curr_POS] = preprocessing(obj,curr_s,curr_h);
            
            %% Data information
            curr_NWindows  = size(curr_P, 1);
            curr_NFreqs    = size(curr_P, 2);
            curr_NChannels = size(curr_P, 3);
            
            
            %% Creating vector labels
            curr_CFeedbackPOS = curr_POS(curr_TYP == 781); % continuos feedback
            curr_CFeedbackDUR = curr_DUR(curr_TYP == 781); % continuos feedback

            curr_CuePOS = curr_POS(curr_TYP == 771 | curr_TYP == 773); % both hands | both feet
            curr_CueDUR = curr_DUR(curr_TYP == 771 | curr_TYP == 773); % both hands | both feet
            curr_CueTYP = curr_TYP(curr_TYP == 771 | curr_TYP == 773); % both hands | both feet

            curr_FixPOS = curr_POS(curr_TYP == 786); % Fixation cross
            curr_FixDUR = curr_DUR(curr_TYP == 786); % Fixation cross
            curr_FixTYP = curr_TYP(curr_TYP == 786); % fixation cross
            
            curr_NumTrials = length(curr_CFeedbackPOS);
            
            % We consider the intersting period from Cue apperance to end of continuous feedback
            curr_Ck = zeros(curr_NWindows, 1);
            curr_Tk = zeros(curr_NWindows, 1);
            curr_TrialStart = nan(curr_NumTrials, 1);
            curr_TrialStop  = nan(curr_NumTrials, 1);
            curr_FixStart = nan(curr_NumTrials, 1);
            curr_FixStop  = nan(curr_NumTrials, 1);
            for trId = 1:curr_NumTrials
                cstart = curr_CuePOS(trId);
                cstop  = curr_CFeedbackPOS(trId) + curr_CFeedbackDUR(trId) - 1;
                curr_Ck(cstart:cstop) = curr_CueTYP(trId);
                curr_Tk(cstart:cstop) = trId;

                curr_TrialStart(trId) = cstart;
                curr_TrialStop(trId)  = cstop;
                curr_FixStart(trId)   = curr_FixPOS(trId);
                curr_FixStop(trId)    = curr_FixPOS(trId) + curr_FixDUR(trId) - 1;
            end
            
            %% Trial extraction

            % Extracting data for each trial (be careful that length might be different for few sample)
            curr_MinTrialDur = min(curr_TrialStop - curr_TrialStart);
            curr_TrialData   = nan(curr_MinTrialDur, curr_NFreqs, curr_NChannels, curr_NumTrials);
            curr_tCk = zeros(curr_NumTrials, 1);
            for trId = 1:curr_NumTrials
                cstart = curr_TrialStart(trId);
                cstop  = cstart + curr_MinTrialDur - 1;
                curr_TrialData(:, :, :, trId)   = curr_P(cstart:cstop, :, :);
                curr_tCk(trId) = unique(curr_Ck(cstart:cstop));
            end
            
            %% Baseline extraction (from fixation)
            curr_MinFixDur = min(curr_FixStop - curr_FixStart);
            curr_FixData   = nan(curr_MinFixDur, curr_NFreqs, curr_NChannels, curr_NumTrials);

            for trId = 1:curr_NumTrials
                cstart = curr_FixStart(trId);
                cstop  = cstart + curr_MinFixDur - 1;
                curr_FixData(:, :, :, trId)   = curr_P(cstart:cstop, :, :);
            end
            
            %% ERD/ERS
            % Average and replicate the value of the baseline
            curr_Baseline = repmat(mean(curr_FixData), [size(curr_TrialData, 1) 1 1 1]);
            curr_ERD = log(curr_TrialData./ curr_Baseline);
            

        end
        
        %% Sessions detector by path
        function obj = listSessions(obj)
            sessions = dir(obj.datapath);
            dfolders = sessions([sessions(:).isdir]);
            dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
            for i = 1 : length(dfolders)
              obj.sessionsNames{i} = dfolders(i).name;
              obj.sessionsPaths{i} = fullfile(obj.datapath,obj.sessionsNames{i});
            end
        end
        
        %% Sessions loader by path
        function obj = loadSessions(obj)
            % First of all load Biosig
            obj.loadBiosig();
            
            for i = 1 : length(obj.sessionsPaths)
                fprintf("Loading session #%d %s \n",i,obj.sessionsPaths{i});
                obj = obj.loadSession(obj.sessionsPaths{i},i);
            end
        end
        
        %% Single session loader by path and number
        function obj = loadSession(obj,path,n)
            gdfPath = fullfile(path,'*.gdf');
            gdfs = dir(gdfPath);
            gdfs = gdfs(~ismember({gdfs(:).name},{'.','..'}));
            
            obj.sessionsData{n}.filenames = [];
            
            for i = 1 : length(gdfs)
                obj.sessionsData{n}.filenames{i} = fullfile(path,gdfs(i).name);
                %fprintf("Loading file #%d %s\n",i,obj.sessionsData{n}.filenames{i});
            end
           
            
            nfiles = length(obj.sessionsData{n}.filenames);
            
            %main merged data
            obj.sessionsData{n}.s   = [];
            obj.sessionsData{n}.TYP = [];
            obj.sessionsData{n}.DUR = [];
            obj.sessionsData{n}.POS = [];
            obj.sessionsData{n}.Rk  = [];
            obj.sessionsData{n}.Mk  = [];
            obj.sessionsData{n}.RkP  = [];
            obj.sessionsData{n}.MkP  = [];
            obj.sessionsData{n}.P     = [];
            obj.sessionsData{n}.freqs = [];
            obj.sessionsData{n}.ERD = [];
            
            obj.allRuns{n}          = 0;
            
            %offline session data merged
            obj.sessionsDataOffline{n}.s   = [];
            obj.sessionsDataOffline{n}.TYP = [];
            obj.sessionsDataOffline{n}.DUR = [];
            obj.sessionsDataOffline{n}.POS = [];
            obj.sessionsDataOffline{n}.Rk  = [];
            obj.sessionsDataOffline{n}.Mk  = [];
            obj.sessionsDataOffline{n}.RkP  = [];
            obj.sessionsDataOffline{n}.MkP  = [];
            obj.sessionsDataOffline{n}.P     = [];
            obj.sessionsDataOffline{n}.freqs = [];
            obj.sessionsDataOffline{n}.ERD = [];
            
            obj.offlineRuns{n}             = 0;
            
            %online data merged
            obj.sessionsDataOnline{n}.s   = [];
            obj.sessionsDataOnline{n}.TYP = [];
            obj.sessionsDataOnline{n}.DUR = [];
            obj.sessionsDataOnline{n}.POS = [];
            obj.sessionsDataOnline{n}.Rk  = [];
            obj.sessionsDataOnline{n}.Mk  = [];
            obj.sessionsDataOnline{n}.RkP  = [];
            obj.sessionsDataOnline{n}.MkP  = [];
            obj.sessionsDataOnline{n}.P     = [];
            obj.sessionsDataOnline{n}.freqs = [];
            obj.sessionsDataOnline{n}.ERD = [];
            obj.onlineRuns{n}             = 0;
             
            for fId=1:nfiles
                % Loading SGD file
                [curr_s, curr_h] = sload(obj.sessionsData{n}.filenames{fId});
                
                % Do preprocessing on single file
                %[curr_P,curr_freqs,curr_TYP,curr_DUR,curr_POS] = obj.preprocessing(curr_s,curr_h);
                
                
                cRk = fId*ones(size(curr_s,1),1);
               
                
                % Chain operation for Offline and Online sets
                if( contains(obj.sessionsData{n}.filenames{fId},'offline') == true)
                    
                    % Do preprocessing on single file
                    [curr_P,curr_freqs,curr_ERD,curr_TYP,curr_DUR,curr_POS] = obj.ERDS(curr_s,curr_h);
                
                    cRkP = fId*ones(size(curr_P,1),1);
                    
                    % Offline session
                    obj.offlineRuns{n} = obj.offlineRuns{n} + 1;
                    cMk = obj.modalityId(1)*ones(size(curr_s,1),1);
                    cMkP = obj.modalityId(1)*ones(size(curr_P,1),1);
                    
                    obj.sessionsDataOffline{n}.TYP = cat(1,obj.sessionsDataOffline{n}.TYP,curr_TYP);
                    obj.sessionsDataOffline{n}.DUR = cat(1,obj.sessionsDataOffline{n}.DUR,curr_DUR);
                    obj.sessionsDataOffline{n}.POS = cat(1,obj.sessionsDataOffline{n}.POS,curr_POS + size(obj.sessionsDataOffline{n}.P,1));
                    obj.sessionsDataOffline{n}.Rk  = cat(1,obj.sessionsDataOffline{n}.Rk,cRk);                  
                    obj.sessionsDataOffline{n}.Mk  = cat(1,obj.sessionsDataOffline{n}.Mk,cMk);
                    obj.sessionsDataOffline{n}.RkP  = cat(1,obj.sessionsDataOffline{n}.RkP,cRkP);                  
                    obj.sessionsDataOffline{n}.MkP  = cat(1,obj.sessionsDataOffline{n}.MkP,cMkP);
                    obj.sessionsDataOffline{n}.s   = cat(1,obj.sessionsDataOffline{n}.s,curr_s);
                    
                                     
                    
                    obj.sessionsDataOffline{n}.P     = cat(1,obj.sessionsDataOffline{n}.P,curr_P);
                    obj.sessionsDataOffline{n}.freqs = cat(1,obj.sessionsDataOffline{n}.freqs,curr_freqs);
                    obj.sessionsDataOffline{n}.ERD{obj.offlineRuns{n}} = curr_ERD;
                    
                    obj.sessionsDataOffline{n}.SampleRate = curr_h.SampleRate;
                    
                    obj.allRuns{n} = obj.allRuns{n} + 1;
                    obj.sessionsData{n}.ERD{obj.allRuns{n}} = curr_ERD; 
                    
                
                elseif( contains(obj.sessionsData{n}.filenames{fId},'online') == true)
                    
                    obj.allRuns{n} = obj.allRuns{n} + 1;
                    
                    % Do preprocessing on single file
                    [curr_P,curr_freqs,curr_TYP,curr_DUR,curr_POS] = obj.preprocessing(curr_s,curr_h);
                
                    cRkP = fId*ones(size(curr_P,1),1);
                    
                    % Online session
                    obj.onlineRuns{n} = obj.onlineRuns{n} + 1;
                    cMk = obj.modalityId(2)*ones(size(curr_s,1),1);
                    cMkP = obj.modalityId(2)*ones(size(curr_P,1),1);
                    
                    obj.sessionsDataOnline{n}.TYP = cat(1,obj.sessionsDataOnline{n}.TYP,curr_TYP);
                    obj.sessionsDataOnline{n}.DUR = cat(1,obj.sessionsDataOnline{n}.DUR,curr_DUR);
                    obj.sessionsDataOnline{n}.POS = cat(1,obj.sessionsDataOnline{n}.POS,curr_POS + size(obj.sessionsDataOnline{n}.P,1));
                    obj.sessionsDataOnline{n}.Rk = cat(1,obj.sessionsDataOnline{n}.Rk,cRk);
                    obj.sessionsDataOnline{n}.Mk = cat(1,obj.sessionsDataOnline{n}.Mk,cMk);
                    obj.sessionsDataOnline{n}.RkP  = cat(1,obj.sessionsDataOnline{n}.RkP,cRkP);                  
                    obj.sessionsDataOnline{n}.MkP  = cat(1,obj.sessionsDataOnline{n}.MkP,cMkP);
                    obj.sessionsDataOnline{n}.s = cat(1,obj.sessionsDataOnline{n}.s,curr_s);
                    
                    obj.sessionsDataOnline{n}.P     = cat(1,obj.sessionsDataOnline{n}.P,curr_P);
                    obj.sessionsDataOnline{n}.freqs = cat(1,obj.sessionsDataOnline{n}.freqs,curr_freqs);
                    %obj.sessionsDataOnline{n}.ERD{obj.onlineRuns{n}} = curr_ERD;
                    
                    obj.sessionsDataOnline{n}.SampleRate = curr_h.SampleRate;
                    
                    
                        
                else
                    %no online - offline data
                    error(['Unknown modality for run: ' obj.sessionsData{n}.filenames{fId}]);
                
                end % online - offline split
                
                %obj.allRuns{n} = obj.allRuns{n} + 1;
                % Chain operation for both type of sessions 
                obj.sessionsData{n}.TYP = cat(1,obj.sessionsData{n}.TYP,curr_TYP);
                obj.sessionsData{n}.DUR = cat(1,obj.sessionsData{n}.DUR,curr_DUR);
                obj.sessionsData{n}.POS = cat(1,obj.sessionsData{n}.POS,curr_POS + size(obj.sessionsData{n}.P,1));              
                obj.sessionsData{n}.Rk  = cat(1,obj.sessionsData{n}.Rk,cRk); 
                obj.sessionsData{n}.RkP  = cat(1,obj.sessionsData{n}.RkP,cRkP);
                obj.sessionsData{n}.s   = cat(1,obj.sessionsData{n}.s,curr_s);
                obj.sessionsData{n}.Mk    = cat(1,obj.sessionsData{n}.Mk,cMk);
                obj.sessionsData{n}.MkP    = cat(1,obj.sessionsData{n}.MkP,cMkP);
                
                obj.sessionsData{n}.P     = cat(1,obj.sessionsData{n}.P,curr_P);
                obj.sessionsData{n}.freqs = cat(1,obj.sessionsData{n}.freqs,curr_freqs);
                %obj.sessionsData{n}.ERD{obj.allRuns{n}} = curr_ERD; 
                obj.sessionsData{n}.SampleRate = curr_h.SampleRate;
                
                
                
            end % session files iteration
            
        end % load session function
        
        
        %% Getter for sessions names
        function sessionsNames = getSessionsNames(obj)
            sessionsNames = obj.sessionsNames;
        end
        
        %% Getter for sessions paths
        function sessionsPaths = getSessionsPaths(obj)
            sessionsPaths = obj.sessionsPaths;
        end
        
        %% Getter for sessions data
        function sessionsData = getSessionsData(obj)
            sessionsData = obj.sessionsData;
        end
        
        %% Getter for offline sessions data
        function offlineSessionsData = getSessionsDataOffline(obj)
            offlineSessionsData = obj.sessionsDataOffline;
        end
        
        %% Getter for online sessions data
        function onlineSessionsData = getSessionsDataOnline(obj)
            onlineSessionsData = obj.sessionsDataOnline;
        end
        
        %% Getter for session data by name
        function sessionData = getSessionByName(obj,name)
            for i=1:length(obj.sessionsNames)
                if name == obj.sessionsNames{i}
                    sessionData = obj.sessionsData{i};
                end
            end
        end
        
        %% Getter for online session data by name
        function sessionDataOnline = getSessionOnlineByName(obj,name)
            for i=1:length(obj.sessionsNames)
                if name == obj.sessionsNames{i}
                    sessionDataOnline = obj.sessionsDataOnline{i};
                end
            end
        end
        
        %% Getter for offline session data by name
        function sessionDataOffline = getSessionOfflineByName(obj,name)
            for i=1:length(obj.sessionsNames)
                if name == obj.sessionsNames{i}
                    sessionDataOffline = obj.sessionsDataOffline{i};
                end
            end
        end
        
        %% Getter for session data by ID
        function sessionData = getSessionById(obj,i)
            sessionData = obj.sessionsData{i};
        end
        
        %% Getter for online session data by ID
        function sessionDataOnline = getSessionOnlineById(obj,i)
            sessionDataOnline = obj.sessionsDataOnline{i};
        end
        
        %% Getter for online session data by ID
        function sessionDataOffline = getSessionOfflineById(obj,i)
            sessionDataOffline = obj.sessionsDataOffline{i};
        end
        
        
    end % methods
    
    %% Static methods
    methods (Static)
        
        %% Biosig loader method
        function loadBiosig()
            %Loading BIOSig Util 3.7.2
            fprintf("Loading Biosig v3.7.2...\n");
            addpath('./Util/biosig4octmat-3.7.2/biosig/t200_FileAccess','./Util/biosig4octmat-3.7.2/biosig/t210_Events','./Util/biosig4octmat-3.7.2/biosig/t250_ArtifactPreProcessingQualityControl','./Util/biosig4octmat-3.7.2/biosig/t300_FeatureExtraction','./Util/biosig4octmat-3.7.2/biosig/t310_ERDSMaps','./Util/biosig4octmat-3.7.2/biosig/t320_Nirs','./Util/biosig4octmat-3.7.2/biosig/t330_StimFit','./Util/biosig4octmat-3.7.2/biosig/t400_Classification','./Util/biosig4octmat-3.7.2/biosig/t450_MultipleTestStatistic','./Util/biosig4octmat-3.7.2/biosig/t490_EvaluationCriteria','./Util/biosig4octmat-3.7.2/biosig/t500_Visualization','./Util/biosig4octmat-3.7.2/biosig/t501_VisualizeCoupling');
        end
    
    end % methods (Static)
end % DataLoader