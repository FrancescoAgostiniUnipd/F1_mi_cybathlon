%% DATA CLASSIFIER FUNCTION CLASS 
%{
Description: Collection of utils and functions used for train 
             and manipulate classifier from the output of data 
             processing

Authors: Agostini Francesco (francesco.agostini.5@studenti.unipd.it)

Version: 1.0

%}

classdef DataClassifier
    %DATACLASSIFIER class
    
    properties
        processor
        % parameters for classifier train
        SelChans
        SelChansId
        SelFreqs
        SelFreqsId
        NumSelFeatures
        LabelIdx
        LabelIdxOnline
        LabelIdxTrain
        LabelIdxTest
        % Single session data for single model
        F
        FOnline
        % Model members
        Models
        
        % Single session data for overall model
        CkTest
        CkTrain
        FTest
        FTrain
        NumWinsTotal
        fullfreqsTotal
        FTotal
        PTotal
        UTotal
        % Model members
        Model
        
        % Accuracy Member
        SSAccTrain
        ppTrain
        GkTrain
        SSClAccTrain
        
        SSAccTest
        ppTest
        ippTest
        GkTest
        SSClAccTest
        
        % Visualization Class
        Presenter
        
    end
    
    methods
        function obj = DataClassifier(p,c,f,pres)%,po_SelChans,po_ChannelLb,po_SelFreqs,po_freqs,po_U,po_Ck,po_NumWins)
            %DATACLASSIFIER Construct an instance of this class (l,p)
            %instances of dataloader and dataprocessor
            
            
            obj.processor       = p;
            %Parameter for classifier training
            obj.SelChans        = c;
            obj.SelFreqs        = f;
            obj.NumSelFeatures  = length(c);
            
            obj.SelChansId      = [];
            obj.SelFreqsId      = [];
            
            obj.F               = [];
            obj.LabelIdx        = [];
            obj.FOnline               = [];
            obj.LabelIdxOnline        = [];
            
            
            obj.CkTest          = [];
            obj.CkTrain         = [];
            obj.FTest           = [];
            obj.FTrain          = [];
            
            obj.NumWinsTotal    = [];
            obj.FTotal          = [];
            obj.UTotal          = [];
            obj.PTotal          = [];
            obj.fullfreqsTotal  = [];
            
            % Accuracy member initialization
            obj.SSAccTrain               = []; 
            obj.ppTrain                  = [];
            obj.GkTrain                  = [];
            obj.SSClAccTrain             = [];
            
            obj.SSAccTest               = []; 
            obj.ppTest                  = [];
            obj.ippTest                 = [];
            obj.GkTest                  = [];
            obj.SSClAccTest             = [];
            
            obj.Presenter       = pres;
            
            obj = obj.loadFromProcessor();
        end
        
        function obj = loadFromProcessor(obj)
            obj = obj.datasetCreator();
            obj = obj.createModel();
            obj = obj.computeTrainsetAccuracy();
            obj = obj.computeTestsetAccuracy();
            obj = obj.computeEvidenceAccumulation();
            obj = obj.computePerformance();
            obj = obj.saveClassifier();
        end
        
        function obj = datasetCreator(obj)
            % Session iterator
            for sId=1:obj.processor.loader.nsessions    
                % Check presence of offline data for trainset creation
                if (obj.processor.loader.offlineRuns{sId} > 0)
                    % IF EXISTS ONLINE RUNS CREATE TRAIN SET
                    [~, SelChansId] = ismember(obj.SelChans, obj.processor.loader.channelLb);
                    [~, SelFreqsId] = ismember(obj.SelFreqs, obj.processor.loader.sessionsDataOffline{sId}.freqs);

                    
                    obj.F{sId} = nan(obj.processor.NumWins{sId}, obj.NumSelFeatures);
                    % Iterate selected features
                    for ftId = 1:obj.NumSelFeatures
                        cfrq  = SelFreqsId(ftId);
                        cchan = SelChansId(ftId);
                        obj.F{sId}(:, ftId) = obj.processor.U{sId}(:, cfrq, cchan);
                    end
                    obj.SelChansId{sId} = SelChansId;
                    obj.SelFreqsId{sId} = SelFreqsId;
                    % Create dataset Labels
                    obj.LabelIdx{sId} = obj.processor.Ck{sId} == 771 | obj.processor.Ck{sId} == 773 | obj.processor.Ck{sId} == 783;
                    % Complete parameter
                    obj.FTrain = cat(1,obj.FTrain,obj.F{sId});
                    obj.CkTrain = cat(1,obj.CkTrain,obj.processor.Ck{sId});
                    
                end  
                
                if (obj.processor.loader.onlineRuns{sId} > 0)
                    % IF EXISTS ONLINE RUNS CREATE TEST SET
                    [~, SelChansId] = ismember(obj.SelChans, obj.processor.loader.channelLb);
                    [~, SelFreqsId] = ismember(obj.SelFreqs, obj.processor.loader.sessionsDataOnline{sId}.freqs);
                    
                    obj.FOnline{sId} = nan(obj.processor.NumWinsOnline{sId}, obj.NumSelFeatures);
                    
                    for ftId = 1:obj.NumSelFeatures
                        cfrq  = SelFreqsId(ftId);
                        cchan = SelChansId(ftId);
                        obj.FOnline{sId}(:, ftId) = obj.processor.UOnline{sId}(:, cfrq, cchan);
                    end
                    
                    obj.LabelIdxOnline{sId} = obj.processor.CkOnline{sId} == 771 | obj.processor.CkOnline{sId} == 773 | obj.processor.CkOnline{sId} == 783;
                    
                    obj.FTest = cat(1,obj.FTest,obj.FOnline{sId});
                    obj.CkTest = cat(1,obj.CkTest,obj.processor.CkOnline{sId});
                end
            end      
        end
        
        function obj = createModel(obj)
            % Session iterator
            for sId=1:obj.processor.loader.nsessions 
                % Check presence of offline data for trainset creation
                if (obj.processor.loader.offlineRuns{sId} > 0)
                    obj.Models{sId} = fitcdiscr(obj.F{sId}(obj.LabelIdx{sId}, :), obj.processor.Ck{sId}(obj.LabelIdx{sId}), 'DiscrimType','quadratic');
                    obj.Presenter.PresentClassifier(obj.processor.loader.sessionsNames{sId},obj.F{sId},obj.processor.Ck{sId},obj.Models{sId},obj.LabelIdx{sId},obj.SelChans,obj.SelFreqs);
                end
            end
            
            obj.LabelIdxTrain = obj.CkTrain == 771 | obj.CkTrain == 773;
            obj.Model = fitcdiscr(obj.FTrain(obj.LabelIdxTrain, :), obj.CkTrain(obj.LabelIdxTrain), 'DiscrimType','quadratic');
            obj.Presenter.PresentClassifier('Global',obj.FTrain,obj.CkTrain,obj.Model,obj.LabelIdxTrain,obj.SelChans,obj.SelFreqs);
        end
        
        function obj = computeTrainsetAccuracy(obj)
            
            for sId=1:obj.processor.loader.nsessions 
                if (obj.processor.loader.offlineRuns{sId} > 0)
                    
                    fprintf("Computing accuracy model #%d\n",sId);
                    NumClasses = length(obj.processor.loader.classId);
                    [obj.GkTrain{sId}, obj.ppTrain{sId}] = predict(obj.Models{sId}, obj.F{sId});

                    obj.SSAccTrain{sId} = 100*sum(obj.GkTrain{sId}(obj.LabelIdx{sId}) == obj.processor.Ck{sId}(obj.LabelIdx{sId}))./length(obj.GkTrain{sId}(obj.LabelIdx{sId}));

                    obj.SSClAccTrain{sId} = nan(NumClasses, 1);
                    for cId = 1:NumClasses
                        cindex = obj.processor.Ck{sId} == obj.processor.loader.classId(cId);
                        obj.SSClAccTrain{sId}(cId) = 100*sum(obj.GkTrain{sId}(cindex) == obj.processor.Ck{sId}(cindex))./length(obj.GkTrain{sId}(cindex));
                    end

                end
                
            end
        end
        
        function obj = computeTestsetAccuracy(obj)
            lastValidModel = 0;
            for sId=1:obj.processor.loader.nsessions 
                if (obj.processor.loader.onlineRuns{sId} > 0) % else no model to test
                    NumClasses = length(obj.processor.loader.classId);
                    if (obj.processor.loader.offlineRuns{sId} > 0)
                        lastValidModel = sId;
                        fprintf("Testing accuracy model #%d\n",sId);      
                        [obj.GkTest{sId}, obj.ppTest{sId}] = predict(obj.Models{sId}, obj.FOnline{sId});
                    else
                        fprintf("Testing accuracy with set #%d, using last known model #%d\n",sId,lastValidModel);
                        [obj.GkTest{sId}, obj.ppTest{sId}] = predict(obj.Models{lastValidModel}, obj.FOnline{sId});
                    end
                    
                    obj.SSAccTest{sId} = 100*sum(obj.GkTest{sId}(obj.LabelIdxOnline{sId}) == obj.processor.CkOnline{sId}(obj.LabelIdxOnline{sId}))./length(obj.GkTest{sId}(obj.LabelIdxOnline{sId}));
                    
                    obj.SSClAccTest{sId} = nan(NumClasses, 1);
                    for cId = 1:NumClasses
                        cindex = obj.processor.CkOnline{sId} == obj.processor.loader.classId(cId);
                        obj.SSClAccTest{sId}(cId) = 100*sum(obj.GkTest{sId}(cindex) == obj.processor.CkOnline{sId}(cindex))./length(obj.GkTest{sId}(cindex));
                    end    
                end
            end
        end
        
        function obj = computeEvidenceAccumulation(obj)
            for sId=1:obj.processor.loader.nsessions 
                if (obj.processor.loader.onlineRuns{sId} > 0)
                    fprintf("Computing Evidence accumulation on model #%d\n",sId);
                    NumClasses = length(obj.processor.loader.classId);
                    TrialStart = obj.processor.loader.sessionsDataOnline{sId}.POS(obj.processor.loader.sessionsDataOnline{sId}.TYP == 781);
                    NumSamples = size(obj.ppTest{sId}, 1);
                    obj.ippTest{sId} = 0.5*ones(size(obj.ppTest{sId}, 1), 1);
                    alpha = 0.97;

                    for samId = 2:NumSamples

                        curr_pp  = obj.ppTest{sId}(samId, 1);
                        prev_ipp = obj.ippTest{sId}(samId-1);

                        if ismember(samId, TrialStart)
                            obj.ippTest{sId}(samId) = 1./NumClasses;
                        else
                            obj.ippTest{sId}(samId) = prev_ipp.*alpha + curr_pp.*(1-alpha);
                        end
                    end
                    
                    % plot result
                    obj.Presenter.PresentAccumulatedEvidence(obj.processor.loader.sessionsNames{sId},obj.processor.TkOnline{sId},obj.processor.CkOnline{sId},obj.ppTest{sId},obj.ippTest{sId},obj.processor.NumTrialsOnline{sId})
                end
                
            end
        end
        
        function obj = computePerformance(obj)
            Threshold     = 0.7;
            for sId=1:obj.processor.loader.nsessions 
                if (obj.processor.loader.onlineRuns{sId} > 0)
                    fprintf("Computing Performance on session #%d\n",sId);
                    ActualClass = obj.processor.loader.sessionsDataOnline{sId}.TYP(obj.processor.loader.sessionsDataOnline{sId}.TYP == 771 | obj.processor.loader.sessionsDataOnline{sId}.TYP == 773 | obj.processor.loader.sessionsDataOnline{sId}.TYP == 783);
                    Decision = nan(obj.processor.NumTrialsOnline{sId}, 1);

                    for trId = 1:obj.processor.NumTrialsOnline{sId}
                        cstart = obj.processor.CFeedbackPOSOnline{sId}(trId);
                        cstop  = obj.processor.CFeedbackPOSOnline{sId}(trId) + obj.processor.CFeedbackDUROnline{sId}(trId) - 1;
                        cipp = obj.ippTest{sId}(cstart:cstop);

                        endpoint = find( (cipp >= Threshold) | (cipp <= 1 - Threshold), 1, 'first' );

                        if(isempty(endpoint))
                            Decision(trId) = 783;
                            continue;
                        end

                        if(cipp(endpoint) >= Threshold)
                            Decision(trId) = 771;
                        elseif (cipp(endpoint) <= Threshold)
                            Decision(trId) = 773;
                        end
                    end
                    
                    ActiveTrials = ActualClass ~= 783;
                    RestTrials = ActualClass == 783;

                    PerfActive  = 100 * (sum(ActualClass(ActiveTrials) == Decision(ActiveTrials))./sum(ActiveTrials));
                    PerfResting = 100 * (sum(ActualClass(RestTrials) == Decision(RestTrials))./sum(RestTrials));

                    RejTrials = Decision == 783;

                    PerfActive_rej = 100 * (sum(ActualClass(ActiveTrials & ~RejTrials) == Decision(ActiveTrials & ~RejTrials))./sum(ActiveTrials & ~RejTrials));
                
                    fprintf("Performance:  Active=%.2f PerfResting=%.2f PerfActiveRej=%.2f\n",PerfActive,PerfResting,PerfActive_rej);
                end
            end
            
            
        end
        
        
        function obj = saveClassifier(obj)
            for sId=1:obj.processor.loader.nsessions
                if (obj.processor.loader.offlineRuns{sId} > 0) % else no model to test
                    Model = obj.Models{sId};
                    SelChansId = obj.SelChansId{sId};
                    SelFreqsId = obj.SelFreqsId{sId};
                    save(strcat('./Classifiers/',obj.processor.loader.sessionsNames{sId}),'Model','SelChansId','SelFreqsId');
                end
            end
        end
        
        function obj = loadClassifier(obj,filename)
            MDL = load(filename);
            obj.SelChansId = MDL.SelChansId;
            obj.SelFreqsId = MDL.SelFreqsId;
            obj.Model = MDL.Model;
        end
    end
end

