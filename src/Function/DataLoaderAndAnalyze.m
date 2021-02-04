%% DATA LOADER FUNCTION CLASS 
%{
Description: Collection of utils and functions userd for load 
              and chain .gdf EEG description file on workspace

Authors: Agostini Francesco (francesco.agostini.5@studenti.unipd.it)
          Deschaux Ophélie   (opheliecandicemarine.deschaux@studenti.unipd.it)
          Marcon Francesco   (francesco.marcon.2@studenti.unipd.it)

Version: 0.1

%}
%% Class Dataloader definition
classdef DataLoaderAndAnalyze
    %% Class Properties 
    properties
        datapath        % Path where data are stored
        datasample      % Sample frequency
        channelLb       % Label fot EEG channels
        channelId       % Channel ID number
        classId         % ID of classes vector
        classLb         % Classes description vector
        nclasses        % N of classes
        modalityId      % ID of modality
        modalityLb      % Label of modality (online ~ offline)
        sessionsNames   % Sessions names present in data directory
        sessionsPaths   % sessions path correspondent to sessions names
        sessionsData    % sessions data correspondent to sessions names
    end
    
    %% Class constructor and loader method
    methods
        %% Constructor
        function obj = DataLoaderAndAnalyze()
            % Setting up env
            obj.datapath   = '../data/';
            obj.datasample = 512;
            obj.channelLb  = {'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'};
            obj.channelId  = 1:length(obj.channelLb);
            obj.classId    = [786 773 771 783 781 897 898];
            obj.classId    = [773 771];
            obj.classLb    = {'Fixation cross','Both hands','Both feet','Rest','Continuos feedback','Target hint','Target miss'}; % we have to add class labels
            obj.classLb    = {'Both hands','Both feet'}; % we have to add class labels;
            obj.nclasses   = length(obj.classId);
            obj.modalityId = [0 1];
            obj.modalityLb = {'offline','online'};
            obj.sessionsData = [];
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

            obj.sessionsData{n}.s=[];
            obj.sessionsData{n}.TYP=[];
            obj.sessionsData{n}.DUR=[];
            obj.sessionsData{n}.POS=[];
            obj.sessionsData{n}.Rk=[];
            obj.sessionsData{n}.Mk=[];
            obj.sessionsData{n}.freqs=[];
            
            mlength    = 1;
            wlength    = 0.5;
            pshift     = 0.25;                  
            wshift     = 0.0625;  
            selfreqs   = 4:2:96;
            winconv = 'backward'; 
             
            for fId=1:nfiles
                %disp(['Loading files (' num2str(fId) '/' num2str(nfiles)'...']);
                [curr_s, curr_h] = sload(obj.sessionsData{n}.filenames{fId});
                
                
                %% Spatial filters
                disp('[proc] |- Applying CAR and Laplacian');
                curr_s = curr_s(:,1:16);
                load('laplacian16.mat');
                s_lap = curr_s*lap;
               

                %% Spectrogram (PSD)
                disp('[proc] |- Computing spectrogram');
                [P, freqgrid] = proc_spectrogram(s_lap, wlength, wshift, pshift, curr_h.SampleRate, mlength);  

                %% Selecting desired frequencies
                [freqs, idfreqs] = intersect(freqgrid, selfreqs);
                P = P(:, idfreqs, :);
                
                %% Extracting events
                disp('[proc] |- Extract and convert the events');
                events.TYP = curr_h.EVENT.TYP;
                events.POS = proc_pos2win(curr_h.EVENT.POS, wshift*curr_h.SampleRate, winconv, mlength*curr_h.SampleRate);
                events.DUR = floor(curr_h.EVENT.DUR/(wshift*curr_h.SampleRate)) + 1;
                events.conversion = winconv;
                
                obj.sessionsData{n}.TYP = cat(1,obj.sessionsData{n}.TYP,events.TYP);
                obj.sessionsData{n}.DUR = cat(1,obj.sessionsData{n}.DUR,events.DUR);
                obj.sessionsData{n}.POS = cat(1,obj.sessionsData{n}.POS,events.POS + size(obj.sessionsData{n}.s,1));

                cRk=fId*ones(size(P,1),1);
                obj.sessionsData{n}.Rk = cat(1,obj.sessionsData{n}.Rk,cRk);

                if( contains(obj.sessionsData{n}.filenames{fId},'offline') == true)
                    cMk = obj.modalityId(1)*ones(size(P,1),1);
                elseif( contains(obj.sessionsData{n}.filenames{fId},'online') == true)
                    cMk = obj.modalityId(2)*ones(size(P,1),1);    
                else
                    error(['Unknown modality for run: ' obj.sessionsData{n}.filenames{fId}]);
                end
                obj.sessionsData{n}.Mk = cat(1,obj.sessionsData{n}.Mk,cMk);

                obj.sessionsData{n}.s = cat(1,obj.sessionsData{n}.s,P);
                obj.sessionsData{n}.SampleRate = curr_h.SampleRate;
                obj.sessionsData{n}.freqs = cat(1,obj.sessionsData{n}.freqs,freqs);

            end
            
        end
        
        
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
        
        %% Getter for session data by name
        function sessionData = getSessionByName(obj,name)
            for i=1:length(obj.sessionsNames)
                if name == obj.sessionsNames{i}
                    sessionData = obj.sessionsData{i};
                end
            end
        end
        
        %% Getter for session data by ID
        function sessionData = getSessionById(obj,i)
            sessionData = obj.sessionsData{i};
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