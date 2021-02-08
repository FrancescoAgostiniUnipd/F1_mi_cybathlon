%% DATA PROCESSING FUNCTION CLASS 
%{
Description: Collection of utils and functions used for process 
             loaded data on workspace

Authors: Agostini Francesco (francesco.agostini.5@studenti.unipd.it)
          Deschaux Ophélie   (opheliecandicemarine.deschaux@studenti.unipd.it)
          Marcon Francesco   (francesco.marcon.2@studenti.unipd.it)

Version: 0.1

%}
%% Class DataProcessing definition
classdef DataProcessing
    
    %% Class Properties 
    properties
        sessions                     % Sessions loaded by dataloader
        sessions_offline             % Sessions Offline
        sessions_online              % Sessions Online
        nsessions                    % Number of sessions
              
    end
    
    methods
        %% Constructor
        function obj = DataProcessing(dataloader)
            % Sessions data loading
            obj.nsessions           = dataloader.getSessionsNumber();
            
            obj.sessions            = dataloader.getSessions();
            obj.sessions_offline    = dataloader.getSessionsOffline();
            obj.sessions_online     = dataloader.getSessionsOnline();
            

        end
        
        %% Sessions Iterator
        function obj = sessionsProcessing()
           for i=1:obj.nsessions
                disp('Processing session n. ',i);
                obj.processSession(i);
            end 
        end
        
        %% Single session processing
        function obj = processSession(index)
            disp('Session Processing...');
            
        end
        
        
        %% Trial Extraction
        function obj = extractTrial(obj,index)
            
        end
        
        
        %% Baseline extraction (from fixation)
        
        
        %% ERD/ERS
        
        %% Labeling the data
        
        %% Fisher Score
        
        
    
        
        
    end %method
end % class