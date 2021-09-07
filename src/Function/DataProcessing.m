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
        % Attributes for Baseline and ERD
        Baseline
        ERD
        % Attributes for trial
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
        MinFixDur
        FixData
        MinTrialDur
        TrialData
        tCk
        SelFreqs
        U
        NumWins
        NumFreqs
        NumChans
        Rk
        Runs
        NumRuns
        freqs
        idfreqs
        NumClasses
        FisherScore
        FS2
        
        
        NWindowsOnline  
        NFreqsOnline   
        NChannelsOnline 
        CFeedbackPOSOnline
        CFeedbackDUROnline
        CuePOSOnline
        CueDUROnline
        CueTYPOnline
        FixPOSOnline
        FixDUROnline
        FixTYPOnline
        NumTrialsOnline
        MinFixDurOnline
        FixDataOnline       
        CkOnline
        TkOnline
        TrialStartOnline
        TrialStopOnline
        FixStartOnline
        FixStopOnline
        MinTrialDurOnline
        TrialDataOnline
        tCkOnline
        UOnline
        NumWinsOnline
        NumFreqsOnline
        NumChansOnline
        RkOnline
        RunsOnline
        NumRunsOnline
        freqsOnline
        idfreqsOnline
        FisherScoreOnline
        FS2Online

        
        % presenter instance
        Presenter
        ForceP
    end
    
    methods
        %% Constructor
        function obj = DataProcessing(dataloader,f,pres)%,forcep)
            %obj.ForceP = forcep;
            
            % Sessions data loading and param
            obj.loader                  = dataloader;
            obj.SelFreqs                = f;
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
            obj.MinFixDur               = [];
            obj.FixData                 = [];
            obj.Baseline                = [];         
            obj.ERD                     = [];
            obj.U                       = [];
            obj.NumWins                 = [];
            obj.NumFreqs                = [];
            obj.NumChans                = [];
            obj.NumClasses              = [];
            obj.Rk                      = [];
            obj.Runs                    = [];
            obj.NumRuns                 = [];
            obj.freqs                   = [];
            obj.idfreqs                 = [];
            obj.FisherScore             = []; 
            obj.FS2                     = [];
            
            
            obj.NWindowsOnline                = [];
            obj.NFreqsOnline                  = [];
            obj.NChannelsOnline               = [];
            obj.CFeedbackPOSOnline            = [];
            obj.CFeedbackDUROnline            = [];
            obj.CuePOSOnline                  = [];
            obj.CueDUROnline                  = [];
            obj.CueTYPOnline                  = [];
            obj.FixPOSOnline                  = [];
            obj.FixDUROnline                  = [];
            obj.FixTYPOnline                  = [];
            obj.CkOnline                      = [];
            obj.TkOnline                      = [];
            obj.NumTrialsOnline               = [];
            obj.TrialStartOnline              = [];
            obj.TrialStopOnline               = [];
            obj.FixStartOnline                = [];
            obj.FixStopOnline                 = [];
            obj.MinTrialDurOnline             = [];
            obj.TrialDataOnline               = [];
            obj.tCkOnline                     = [];
            % Baseline data
    
            obj.UOnline                 = [];
            obj.NumWinsOnline           = [];
            obj.NumFreqsOnline          = [];
            obj.NumChansOnline          = [];
            obj.RkOnline                = [];
            obj.RunsOnline              = [];
            obj.NumRunsOnline           = [];
            obj.idfreqsOnline           = [];
            
            obj.FisherScoreOnline             = []; 
            obj.FS2Online                     = [];
            
            obj.Presenter = pres;
            
            % Process session from dataloader
            obj = obj.sessionsProcessing();
            
        end
        
        %% Sessions Iterator
        function obj = sessionsProcessing(obj)
           for i=1:obj.loader.nsessions
                fprintf("Processing session #%d %s \n",i,obj.loader.sessionsPaths{i});
                obj = obj.processSession(i);
            end 
        end
        
        %% Single session processing
        function obj = processSession(obj,index)
            obj = obj.baseVector(index);
            obj = obj.extractTrial(index);
            obj = obj.baselineExtraction(index);
            obj = obj.ErsErdComputing(index);
            obj = obj.LogData(index);
            obj = obj.LabelingData(index);
            obj = obj.computeFisherScore(index);
        end
        
        
        %% Base Vector
        function obj = baseVector(obj,i)
            if (obj.loader.onlineRuns{i} > 0)
                % Base vectors
                obj.NWindowsOnline{i}  = size(obj.loader.sessionsDataOnline{i}.P, 1);
                obj.NFreqsOnline{i}    = size(obj.loader.sessionsDataOnline{i}.P, 2);
                obj.NChannelsOnline{i} = size(obj.loader.sessionsDataOnline{i}.P, 3);
                obj.CFeedbackPOSOnline{i} = obj.loader.sessionsDataOnline{i}.POS(obj.loader.sessionsDataOnline{i}.TYP == 781);
                obj.CFeedbackDUROnline{i} = obj.loader.sessionsDataOnline{i}.DUR(obj.loader.sessionsDataOnline{i}.TYP == 781);
                obj.CuePOSOnline{i} = obj.loader.sessionsDataOnline{i}.POS(obj.loader.sessionsDataOnline{i}.TYP == 771 | obj.loader.sessionsDataOnline{i}.TYP == 773 | obj.loader.sessionsDataOnline{i}.TYP == 783 );
                obj.CueDUROnline{i} = obj.loader.sessionsDataOnline{i}.DUR(obj.loader.sessionsDataOnline{i}.TYP == 771 | obj.loader.sessionsDataOnline{i}.TYP == 773 | obj.loader.sessionsDataOnline{i}.TYP == 783 );
                obj.CueTYPOnline{i} = obj.loader.sessionsDataOnline{i}.TYP(obj.loader.sessionsDataOnline{i}.TYP == 771 | obj.loader.sessionsDataOnline{i}.TYP == 773 | obj.loader.sessionsDataOnline{i}.TYP == 783 );
                obj.FixPOSOnline{i} = obj.loader.sessionsDataOnline{i}.POS(obj.loader.sessionsDataOnline{i}.TYP == 786);
                obj.FixDUROnline{i} = obj.loader.sessionsDataOnline{i}.DUR(obj.loader.sessionsDataOnline{i}.TYP == 786);
                obj.FixTYPOnline{i} = obj.loader.sessionsDataOnline{i}.TYP(obj.loader.sessionsDataOnline{i}.TYP == 786);
                obj.NumTrialsOnline{i} = length(obj.CFeedbackPOSOnline{i});
            end
            
            if (obj.loader.offlineRuns{i} > 0)
                % Base vectors
                obj.NWindows{i}  = size(obj.loader.sessionsDataOffline{i}.P, 1);
                obj.NFreqs{i}    = size(obj.loader.sessionsDataOffline{i}.P, 2);
                obj.NChannels{i} = size(obj.loader.sessionsDataOffline{i}.P, 3);
                obj.CFeedbackPOS{i} = obj.loader.sessionsDataOffline{i}.POS(obj.loader.sessionsDataOffline{i}.TYP == 781);
                obj.CFeedbackDUR{i} = obj.loader.sessionsDataOffline{i}.DUR(obj.loader.sessionsDataOffline{i}.TYP == 781);
                obj.CuePOS{i} = obj.loader.sessionsDataOffline{i}.POS(obj.loader.sessionsDataOffline{i}.TYP == 771 | obj.loader.sessionsDataOffline{i}.TYP == 773 | obj.loader.sessionsDataOffline{i}.TYP == 783 );
                obj.CueDUR{i} = obj.loader.sessionsDataOffline{i}.DUR(obj.loader.sessionsDataOffline{i}.TYP == 771 | obj.loader.sessionsDataOffline{i}.TYP == 773 | obj.loader.sessionsDataOffline{i}.TYP == 783 );
                obj.CueTYP{i} = obj.loader.sessionsDataOffline{i}.TYP(obj.loader.sessionsDataOffline{i}.TYP == 771 | obj.loader.sessionsDataOffline{i}.TYP == 773 | obj.loader.sessionsDataOffline{i}.TYP == 783 );
                obj.FixPOS{i} = obj.loader.sessionsDataOffline{i}.POS(obj.loader.sessionsDataOffline{i}.TYP == 786);
                obj.FixDUR{i} = obj.loader.sessionsDataOffline{i}.DUR(obj.loader.sessionsDataOffline{i}.TYP == 786);
                obj.FixTYP{i} = obj.loader.sessionsDataOffline{i}.TYP(obj.loader.sessionsDataOffline{i}.TYP == 786);
                obj.NumTrials{i} = length(obj.CFeedbackPOS{i});
            end
            
           
        end
        %% Trial Extraction
        function obj = extractTrial(obj,i)

            if (obj.loader.onlineRuns{i} > 0)
                obj.CkOnline{i} = zeros(obj.NWindowsOnline{i}, 1);
                obj.TkOnline{i} = zeros(obj.NWindowsOnline{i}, 1);
                obj.TrialStartOnline{i} = nan(obj.NumTrialsOnline{i}, 1);
                obj.TrialStopOnline{i}  = nan(obj.NumTrialsOnline{i}, 1);
                obj.FixStartOnline{i} = nan(obj.NumTrialsOnline{i}, 1);
                obj.FixStopOnline{i}  = nan(obj.NumTrialsOnline{i}, 1);
                for trId = 1:obj.NumTrialsOnline{i}
                    cstart = obj.CuePOSOnline{i}(trId);
                    cstop  = obj.CFeedbackPOSOnline{i}(trId) + obj.CFeedbackDUROnline{i}(trId) - 1;
                    obj.CkOnline{i}(cstart:cstop) = obj.CueTYPOnline{i}(trId);
                    obj.TkOnline{i}(cstart:cstop) = trId;

                    obj.TrialStartOnline{i}(trId) = cstart;
                    obj.TrialStopOnline{i}(trId)  = cstop;
                    obj.FixStartOnline{i}(trId)   = obj.FixPOSOnline{i}(trId);
                    obj.FixStopOnline{i}(trId)    = obj.FixPOSOnline{i}(trId) + obj.FixDUROnline{i}(trId) - 1;
                end

                % extraction
                obj.MinTrialDurOnline{i} = min(obj.TrialStopOnline{i} - obj.TrialStartOnline{i});
                %fprintf("i = %d | MTD{i} = %f | TStart{i} = %f | TStop{i} = %f \n",i,obj.MinTrialDur{i},obj.TrialStop{i},obj.TrialStart{i});
                obj.TrialDataOnline{i}   = nan(obj.MinTrialDurOnline{i}, obj.NFreqsOnline{i}, obj.NChannelsOnline{i}, obj.NumTrialsOnline{i});
                obj.tCkOnline{i} = zeros(obj.NumTrialsOnline{i}, 1);
                for trId = 1:obj.NumTrialsOnline{i}
                    cstart = obj.TrialStartOnline{i}(trId);
                    cstop  = cstart + obj.MinTrialDurOnline{i} - 1;
                    obj.TrialDataOnline{i}(:, :, :, trId)   = obj.loader.sessionsDataOnline{i}.P(cstart:cstop, :, :);
                    obj.tCkOnline{i}(trId) = unique(obj.CkOnline{i}(cstart:cstop));
                end
            end    
            
            if (obj.loader.offlineRuns{i} > 0)
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
                %fprintf("i = %d | MTD{i} = %f | TStart{i} = %f | TStop{i} = %f \n",i,obj.MinTrialDur{i},obj.TrialStop{i},obj.TrialStart{i});
                obj.TrialData{i}   = nan(obj.MinTrialDur{i}, obj.NFreqs{i}, obj.NChannels{i}, obj.NumTrials{i});
                obj.tCk{i} = zeros(obj.NumTrials{i}, 1);
                for trId = 1:obj.NumTrials{i}
                    cstart = obj.TrialStart{i}(trId);
                    cstop  = cstart + obj.MinTrialDur{i} - 1;
                    obj.TrialData{i}(:, :, :, trId)   = obj.loader.sessionsDataOffline{i}.P(cstart:cstop, :, :);

                    obj.tCk{i}(trId) = unique(obj.Ck{i}(cstart:cstop));
                end
            end
            

        end
        
        
        %% Baseline extraction
        
        function obj = baselineExtraction(obj,i)
            if (obj.loader.onlineRuns{i} > 0)
                obj.MinFixDurOnline{i} = min(obj.FixStopOnline{i} - obj.FixStartOnline{i});
                obj.FixDataOnline{i}   = nan(obj.MinFixDurOnline{i}, obj.NFreqsOnline{i}, obj.NChannelsOnline{i}, obj.NumTrialsOnline{i});
                for trId = 1:obj.NumTrialsOnline{i}
                    cstart = obj.FixStartOnline{i}(trId);
                    cstop  = cstart + obj.MinFixDurOnline{i} - 1;
                    obj.FixDataOnline{i}(:, :, :, trId)   = obj.loader.sessionsDataOnline{i}.P(cstart:cstop, :, :);
                end 
            end
            
            if (obj.loader.offlineRuns{i} > 0)
                obj.MinFixDur{i} = min(obj.FixStop{i} - obj.FixStart{i});
                obj.FixData{i}   = nan(obj.MinFixDur{i}, obj.NFreqs{i}, obj.NChannels{i}, obj.NumTrials{i});
                for trId = 1:obj.NumTrials{i}
                    cstart = obj.FixStart{i}(trId);
                    cstop  = cstart + obj.MinFixDur{i} - 1;
                    obj.FixData{i}(:, :, :, trId)   = obj.loader.sessionsDataOffline{i}.P(cstart:cstop, :, :);
                end 
            end
        end
        
        
        %% ERD/ERS
        
        function obj = ErsErdComputing(obj,i)
            if (obj.loader.offlineRuns{i} > 0)
                obj.Baseline{i} = repmat(mean(obj.FixData{i}), [size(obj.TrialData{i}, 1) 1 1 1]);
                obj.ERD{i} = log(obj.TrialData{i}./ obj.Baseline{i}); 
            end
        end    
        
        %% Labeling the data
        
        function obj = LogData(obj,i)
            
            if (obj.loader.onlineRuns{i} > 0)
                % Online data
                fullFreqsOnline = obj.loader.sessionsDataOnline{i}.freqs;
                [obj.freqsOnline{i}, obj.idfreqsOnline{i}] = intersect(fullFreqsOnline, obj.SelFreqs);
                obj.UOnline{i} = log(obj.loader.sessionsDataOnline{i}.P(:, obj.idfreqsOnline{i}, :));
                obj.NumWinsOnline{i}  = size(obj.UOnline{i}, 1);
                obj.NumFreqsOnline{i} = size(obj.UOnline{i}, 2);
                obj.NumChansOnline{i} = size(obj.UOnline{i}, 3);
                obj.RkOnline{i} = obj.loader.sessionsDataOnline{i}.RkP;
                obj.RunsOnline{i} = unique(obj.RkOnline{i});
                obj.NumRunsOnline{i} = length(obj.RunsOnline{i});
            end
            if (obj.loader.offlineRuns{i} > 0)
                % offline data
                fullFreqs = obj.loader.sessionsDataOffline{i}.freqs;
                [obj.freqs{i}, obj.idfreqs{i}] = intersect(fullFreqs, obj.SelFreqs);

                obj.U{i} = log(obj.loader.sessionsDataOffline{i}.P(:, obj.idfreqs{i}, :));
                %obj.U{i} = log(obj.ForceP(:, obj.idfreqs{i}, :));

                obj.NumWins{i}  = size(obj.U{i}, 1);
                obj.NumFreqs{i} = size(obj.U{i}, 2);
                obj.NumChans{i} = size(obj.U{i}, 3);

                obj.Rk{i} = obj.loader.sessionsDataOffline{i}.RkP;
                obj.Runs{i} = unique(obj.Rk{i});

                obj.NumRuns{i} = length(obj.Runs{i});
            end
            
        end

        function obj = LabelingData(obj,i)
            if (obj.loader.onlineRuns{i} > 0)
                % Online data
                obj.CkOnline{i} = zeros(obj.NumWinsOnline{i}, 1);
                obj.TkOnline{i} = zeros(obj.NumWinsOnline{i}, 1);
                for trId = 1:obj.NumTrialsOnline{i}
                    cstart = obj.CuePOSOnline{i}(trId);
                    cstop  = obj.CFeedbackPOSOnline{i}(trId) + obj.CFeedbackDUROnline{i}(trId) - 1;
                    obj.CkOnline{i}(cstart:cstop) = obj.CueTYPOnline{i}(trId);
                    obj.TkOnline{i}(cstart:cstop) = trId;
                end
            end
            if (obj.loader.offlineRuns{i} > 0)
                % Offline data
                obj.Ck{i} = zeros(obj.NumWins{i}, 1);
                obj.Tk{i} = zeros(obj.NumWins{i}, 1);
                for trId = 1:obj.NumTrials{i}
                    cstart = obj.CuePOS{i}(trId);
                    cstop  = obj.CFeedbackPOS{i}(trId) + obj.CFeedbackDUR{i}(trId) - 1;
                    obj.Ck{i}(cstart:cstop) = obj.CueTYP{i}(trId);
                    obj.Tk{i}(cstart:cstop) = trId;
                end
            end
            
        end
        
        %% Fisher Score
        
        function obj =  computeFisherScore(obj,i)
            
            % TODO ONLINE
            
            
            if (obj.loader.offlineRuns{i} > 0)
                obj.NumClasses{i} = length(obj.loader.classId);

                obj.FisherScore{i} = nan(obj.NumFreqs{i}, obj.NumChans{i}, obj.NumRuns{i});
                obj.FS2{i} = nan(obj.NumFreqs{i}*obj.NumChans{i}, obj.NumRuns{i});
                for rId = 1:obj.NumRuns{i}
                    rindex = obj.Rk{i} == obj.Runs{i}(rId); 

                    cmu    = nan(obj.NumFreqs{i}, obj.NumChans{i}, 2);
                    csigma = nan(obj.NumFreqs{i}, obj.NumChans{i}, 2);

                    for cId = 1:obj.NumClasses{i}
                        cindex = rindex & obj.Ck{i} == obj.loader.classId(cId);
                        cmu(:, :, cId) = squeeze(mean(obj.U{i}(cindex, :, :)));
                        csigma(:, :, cId) = squeeze(std(obj.U{i}(cindex, :, :)));
                    end

                    obj.FisherScore{i}(:, :, rId) = abs(cmu(:, :, 2) - cmu(:, :, 1)) ./ sqrt( ( csigma(:, :, 1).^2 + csigma(:, :, 2).^2 ) );
                end

                % print fisher score debug

                if (obj.loader.offlineRuns{i} == 0)
                    %fprintf("No offline to visualize for session %d\n",i);
                else
                    %obj.Presenter.PresentErdErs(obj.loader.sessionsNames{i}, obj.MinTrialDur{i}, obj.loader.wshift, 0 ,obj.NumClasses{i},obj.ERD{i},obj.loader.classId,obj.freqs{i},obj.loader.channelLb,obj.loader.classLb,obj.tCk{i});
                    obj.Presenter.PresentFisherScore(obj.loader.sessionsNames{i}, obj.NumRuns{i},obj.NumFreqs{i}, obj.NumChans{i}, obj.loader.channelLb, obj.freqs{i}, obj.FisherScore{i})
                end
            end
        end % ComputeFisherScore
        
    
        
        
    end %method
end % class