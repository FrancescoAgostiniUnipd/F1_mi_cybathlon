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
        loader                       % DataLoader refs
        
        NWindows  
        NFreqs   
        NChannels 
        CFeedbackPOS
        CFeedbackDUR
        CuePOS
        CueDUR
        CueTYP
        FixPOS
        FixDUR
        FixTYP
        NumTrials
        Ck
        Tk
        TrialStart
        TrialStop
        FixStart
        FixStop
        MinTrialDur
        TrialData
        tCk

        
    end
    
    methods
        %% Constructor
        function obj = DataProcessing(dataloader)
            % Sessions data loading
            obj.loader              = dataloader;
            % Trial Vectors
            NWindows                = [];
            NFreqs                  = [];
            NChannels               = [];
            CFeedbackPOS            = [];
            CFeedbackDUR            = [];
            CuePOS                  = [];
            CueDUR                  = [];
            CueTYP                  = [];
            FixPOS                  = [];
            FixDUR                  = [];
            FixTYP                  = [];
            NumTrials               = [];
            Ck                      = [];
            Tk                      = [];
            TrialStart              = [];
            TrialStop               = [];
            FixStart                = [];
            FixStop                 = [];
            MinTrialDur             = [];
            TrialData               = [];
            tCk                     = [];

            
            obj.sessionsProcessing();
        end
        
        %% Sessions Iterator
        function obj = sessionsProcessing(obj)
           for i=1:obj.loader.nsessions
                fprintf("Processing session #%d %s \n",i,obj.loader.sessionsPaths{i});
                obj.processSession(i);
            end 
        end
        
        %% Single session processing
        function obj = processSession(obj,index)
            disp('Session Processing...');
            obj.extractTrial(index);
        end
        
        
        %% Trial Extraction
        function obj = extractTrial(obj,i)
            
            % Base vectors
            obj.NWindows{i}  = size(obj.loader.sessionsDataOffline{i}.P, 1);
            obj.NFreqs{i}    = size(obj.loader.sessionsDataOffline{i}.P, 2);
            obj.NChannels{i} = size(obj.loader.sessionsDataOffline{i}.P, 3);
            obj.CFeedbackPOS{i} = obj.loader.sessionsDataOffline{i}.POS(obj.loader.sessionsDataOffline{i}.TYP == 781);
            obj.CFeedbackDUR{i} = obj.loader.sessionsDataOffline{i}.DUR(obj.loader.sessionsDataOffline{i}.TYP == 781);
            obj.CuePOS{i} = obj.loader.sessionsDataOffline{i}.POS(obj.loader.sessionsDataOffline{i}.TYP == 771 | obj.loader.sessionsDataOffline{i}.TYP == 773 );
            obj.CueDUR{i} = obj.loader.sessionsDataOffline{i}.DUR(obj.loader.sessionsDataOffline{i}.TYP == 771 | obj.loader.sessionsDataOffline{i}.TYP == 773 );
            obj.CueTYP{i} = obj.loader.sessionsDataOffline{i}.TYP(obj.loader.sessionsDataOffline{i}.TYP == 771 | obj.loader.sessionsDataOffline{i}.TYP == 773);
            obj.FixPOS{i} = obj.loader.sessionsDataOffline{i}.POS(obj.loader.sessionsDataOffline{i}.TYP == 786);
            obj.FixDUR{i} = obj.loader.sessionsDataOffline{i}.DUR(obj.loader.sessionsDataOffline{i}.TYP == 786);
            obj.FixTYP{i} = obj.loader.sessionsDataOffline{i}.TYP(obj.loader.sessionsDataOffline{i}.TYP == 786);
            obj.NumTrials{i} = length(obj.CFeedbackPOS{i});
            
            % Trial 
            obj.Ck{i} = zeros(obj.NWindows{i}, 1);
            obj.Tk{i} = zeros(obj.NWindows{i}, 1);
            obj.TrialStart{i} = nan(obj.NumTrials{i}, 1);
            obj.TrialStop{i}  = nan(obj.NumTrials{i}, 1);
            obj.FixStart{i} = nan(obj.NumTrials{i}, 1);
            obj.FixStop{i}  = nan(obj.NumTrials{i}, 1);
            for trId = 1:obj.NumTrials{i}
                cstart = obj.CuePOS{i}(trId);
                cstop  = obj.CFeedbackPOS{i}(trId) + obj.CFeedbackDUR{i}(trId) - 1;
                obj.Ck{i}(cstart:cstop) = obj.CueTYP{i}(trId);
                obj.Tk{i}(cstart:cstop) = trId;

                obj.TrialStart{i}(trId) = cstart;
                obj.TrialStop{i}(trId)  = cstop;
                obj.FixStart{i}(trId)   = obj.FixPOS{i}(trId);
                obj.FixStop{i}(trId)    = obj.FixPOS{i}(trId) + obj.FixDUR{i}(trId) - 1;
            end

            % extraction
            obj.MinTrialDur{i} = min(obj.TrialStop{i} - obj.TrialStart{i});
            obj.TrialData{i}   = nan(obj.MinTrialDur{i}, obj.NFreqs{i}, obj.NChannels{i}, obj.NumTrials{i});
            obj.tCk{i} = zeros(obj.NumTrials{i}, 1);
            for trId = 1:obj.NumTrials{i}
                cstart = obj.TrialStart{i}(trId);
                cstop  = cstart + obj.MinTrialDur{i} - 1;
                obj.TrialData{i}(:, :, :, trId)   = obj.loader.sessionsDataOffline{i}.P(cstart:cstop, :, :);

                obj.tCk{i}(trId) = unique(obj.Ck{i}(cstart:cstop));
            end
        end
        
        
        %% Baseline extraction (from fixation)
        
        
        %% ERD/ERS
        
        %% Labeling the data
        
        %% Fisher Score
        
        
    
        
        
    end %method
end % class