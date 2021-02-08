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
        
        MinFixDur
        FixData
        
        Baseline
        ERD

        SelFreqs
    end
    
    methods
        %% Constructor
        function obj = DataProcessing(dataloader)
            % Sessions data loading
            obj.loader              = dataloader;
            % Trial Vectors
            obj.NWindows                = [];
            obj.NFreqs                  = [];
            obj.NChannels               = [];
            obj.CFeedbackPOS            = [];
            obj.CFeedbackDUR            = [];
            obj.CuePOS                  = [];
            obj.CueDUR                  = [];
            obj.CueTYP                  = [];
            obj.FixPOS                  = [];
            obj.FixDUR                  = [];
            obj.FixTYP                  = [];
            obj.NumTrials               = [];
            obj.Ck                      = [];
            obj.Tk                      = [];
            obj.TrialStart              = [];
            obj.TrialStop               = [];
            obj.FixStart                = [];
            obj.FixStop                 = [];
            obj.MinTrialDur             = [];
            obj.TrialData               = [];
            obj.tCk                     = [];
            % Baseline data
            obj.MinFixDur               = [];
            obj.FixData                 = [];
            obj.baseline                = [];
             
            obj.ERD                     = [];
            
            obj.SelFreqs = 4:2:48;
            
            
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
            obj = obj.baseVector(index);
            obj = obj.extractTrial(index);
            obj = obj.baselineExtraction(index);
            obj = obj.ErsErdComputing(index);
        end
        
        
        %% Base Vector
        function obj = baseVector(obj,i)
            
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
           
           
        end
        %% Trial Extraction
        function obj = extractTrial(obj,i)

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
        
        
        %% Baseline extraction
        
        function obj = baselineExtraction(obj,i)
            
            obj.MinFixDur{i} = min(obj.FixStop{i} - obj.FixStart{i});
            obj.FixData{i}   = nan(obj.MinFixDur{i}, obj.NFreqs{i}, obj.NChannels{i}, obj.NumTrials{i});
                   
            for trId = 1:obj.NumTrials{i}
                cstart = obj.FixStart{i}(trId);
                cstop  = cstart + obj.MinFixDur{i} - 1;
                obj.FixData{i}(:, :, :, trId)   = obj.loader.sessionsDataOffline{i}.P(cstart:cstop, :, :);
            end 
        end
        
        
        %% ERD/ERS
        
        function obj = ErsErdComputing(obj,i)
            obj.Baseline{i} = repmat(mean(obj.FixData{i}), [size(obj.TrialData{i}, 1) 1 1 1]);
            obj.ERD{i} = log(obj.TrialData{i}./ obj.Baseline{i}); 
        end    
        
        %% Labeling the data
        
        function obj = LogData(obj,i)
            
            fullFreqs = obj.loader.sessionsDataOffline{i}.freqs;
            [freqs, idfreqs] = intersect(fullFreqs, SelFreqs);

            U = log(sessionOffline1.P);

            NumWins  = size(U, 1);
            NumFreqs = size(U, 2);
            NumChans = size(U, 3);

            Rk = obj.loader.sessionsDataOffline{i}.RkP;
            Runs = unique(Rk);

            NumRuns = length(Runs);
            
        end

        function obj = LabelingData(obj,i)

            % 
            % We consider the intersting period from Cue apperance to end of continuous feedback
            Ck = zeros(NumWins, 1);
            Tk = zeros(NumWins, 1);

            for trId = 1:NumTrials
                cstart = CuePOS(trId);
                cstop  = CFeedbackPOS(trId) + CFeedbackDUR(trId) - 1;
                Ck(cstart:cstop) = CueTYP(trId);
                Tk(cstart:cstop) = trId;
            end
    
        end
        
        %% Fisher Score
        
        function obj =  computeFisherScore(obj,i)

            NumClasses = length(data.classId);

            FisherScore = nan(NumFreqs, NumChans, NumRuns);
            FS2 = nan(NumFreqs*NumChans, NumRuns);
            for rId = 1:NumRuns
                rindex = Rk == Runs(rId); 

                cmu    = nan(NumFreqs, NumChans, 2);
                csigma = nan(NumFreqs, NumChans, 2);

                for cId = 1:NumClasses
                    cindex = rindex & Ck == data.classId(cId);
                    cmu(:, :, cId) = squeeze(mean(U(cindex, :, :)));
                    csigma(:, :, cId) = squeeze(std(U(cindex, :, :)));
                end

                FisherScore(:, :, rId) = abs(cmu(:, :, 2) - cmu(:, :, 1)) ./ sqrt( ( csigma(:, :, 1).^2 + csigma(:, :, 2).^2 ) );
            end
            
        end % ComputeFisherScore
        
    
        
        
    end %method
end % class