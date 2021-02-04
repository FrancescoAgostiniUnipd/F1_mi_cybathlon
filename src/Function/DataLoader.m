%% DATA LOADER FUNCTION CLASS 
%{
Description: Collection of utils and functions userd for load 
              and chain .gdf EEG description file on workspace

Authors: Agostini Francesco (francesco.agostini.5@studenti.unipd.it)
          Deschaux Ophélie   (opheliecandicemarine.deschaux@studenti.unipd.it)
          Marcon Francesco   (francesco.marcon.2@studenti.unipd.it)

Version: 0.3

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
            obj = obj.loadSessions();
            
            sprintf("\nAll %d sessions loaded correctly!\n\n",length(obj.sessionsNames));
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
            obj.allRuns{n}          = 0;
            
            %offline session data merged
            obj.sessionsDataOffline{n}.s   = [];
            obj.sessionsDataOffline{n}.TYP = [];
            obj.sessionsDataOffline{n}.DUR = [];
            obj.sessionsDataOffline{n}.POS = [];
            obj.sessionsDataOffline{n}.Rk  = [];
            obj.sessionsDataOffline{n}.Mk  = [];
            obj.offlineRuns{n}             = 0;
            
            %online data merged
            obj.sessionsDataOnline{n}.s   = [];
            obj.sessionsDataOnline{n}.TYP = [];
            obj.sessionsDataOnline{n}.DUR = [];
            obj.sessionsDataOnline{n}.POS = [];
            obj.sessionsDataOnline{n}.Rk  = [];
            obj.sessionsDataOnline{n}.Mk  = [];
            obj.onlineRuns{n}             = 0;
             
            for fId=1:nfiles
                % Loading SGD file
                [curr_s, curr_h] = sload(obj.sessionsData{n}.filenames{fId});
                
                cRk = fId*ones(size(curr_s,1),1);
                
                % Chain operation for Offline and Online sets
                if( contains(obj.sessionsData{n}.filenames{fId},'offline') == true)
                    
                    % Offline session
                    cMk = obj.modalityId(1)*ones(size(curr_s,1),1);
                    
                    obj.sessionsDataOffline{n}.TYP = cat(1,obj.sessionsDataOffline{n}.TYP,curr_h.EVENT.TYP);
                    obj.sessionsDataOffline{n}.DUR = cat(1,obj.sessionsDataOffline{n}.DUR,curr_h.EVENT.DUR);
                    obj.sessionsDataOffline{n}.POS = cat(1,obj.sessionsDataOffline{n}.POS,curr_h.EVENT.POS + size(obj.sessionsDataOffline{n}.s,1));
                    obj.sessionsDataOffline{n}.Rk  = cat(1,obj.sessionsDataOffline{n}.Rk,cRk);                  
                    obj.sessionsDataOffline{n}.Mk  = cat(1,obj.sessionsDataOffline{n}.Mk,cMk);
                    obj.sessionsDataOffline{n}.s   = cat(1,obj.sessionsDataOffline{n}.s,curr_s);
                    obj.sessionsDataOffline{n}.SampleRate = curr_h.SampleRate;
                    obj.offlineRuns{n} = obj.offlineRuns{n} + 1;
                    
                    
                
                elseif( contains(obj.sessionsData{n}.filenames{fId},'online') == true)
                    
                    % Online session
                    cMk = obj.modalityId(2)*ones(size(curr_s,1),1);
                    
                    obj.sessionsDataOnline{n}.TYP = cat(1,obj.sessionsDataOnline{n}.TYP,curr_h.EVENT.TYP);
                    obj.sessionsDataOnline{n}.DUR = cat(1,obj.sessionsDataOnline{n}.DUR,curr_h.EVENT.DUR);
                    obj.sessionsDataOnline{n}.POS = cat(1,obj.sessionsDataOnline{n}.POS,curr_h.EVENT.POS + size(obj.sessionsDataOnline{n}.s,1));
                    obj.sessionsDataOnline{n}.Rk = cat(1,obj.sessionsDataOnline{n}.Rk,cRk);
                    obj.sessionsDataOnline{n}.Mk = cat(1,obj.sessionsDataOnline{n}.Mk,cMk);
                    obj.sessionsDataOnline{n}.s = cat(1,obj.sessionsDataOnline{n}.s,curr_s);
                    obj.sessionsDataOnline{n}.SampleRate = curr_h.SampleRate;
                    obj.onlineRuns{n} = obj.onlineRuns{n} + 1;
                    
                        
                else
                    %no online - offline data
                    error(['Unknown modality for run: ' obj.sessionsData{n}.filenames{fId}]);
                
                end % online - offline split
                
                % Chain operation for both type of sessions 
                obj.sessionsData{n}.TYP = cat(1,obj.sessionsData{n}.TYP,curr_h.EVENT.TYP);
                obj.sessionsData{n}.DUR = cat(1,obj.sessionsData{n}.DUR,curr_h.EVENT.DUR);
                obj.sessionsData{n}.POS = cat(1,obj.sessionsData{n}.POS,curr_h.EVENT.POS + size(obj.sessionsData{n}.s,1));              
                obj.sessionsData{n}.Rk  = cat(1,obj.sessionsData{n}.Rk,cRk);                            
                obj.sessionsData{n}.s   = cat(1,obj.sessionsData{n}.s,curr_s);
                obj.sessionsData{n}.SampleRate = curr_h.SampleRate;
                obj.allRuns{n} = obj.allRuns{n} + 1;
                obj.sessionsData{n}.Mk  = cat(1,obj.sessionsData{n}.Mk,cMk);
                
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