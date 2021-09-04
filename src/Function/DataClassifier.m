%% DATA CLASSIFIER FUNCTION CLASS 
%{
Description: Collection of utils and functions used for train 
             and manipulate classifier from the output of data 
             processing

Authors: Agostini Francesco (francesco.agostini.5@studenti.unipd.it)
          Deschaux Ophélie   (opheliecandicemarine.deschaux@studenti.unipd.it)
          Marcon Francesco   (francesco.marcon.2@studenti.unipd.it)

Version: 0.1

%}

classdef DataClassifier
    %DATACLASSIFIER class
    
    properties
        processor
        % parameters for classifier train
        SelChans
        SelFreqs
        NumSelFeatures
        LabelIdx
        SelChansId
        % Single session data for single model
        Ck
        NumWins
        fullfreqs
        F
        P
        U
        % Model members
        Models
        
        % Single session data for overall model
        CkTotal
        NumWinsTotal
        fullfreqsTotal
        FTotal
        PTotal
        UTotal
        % Model members
        Model
        
        % Accuracy Member
        SSAcc
        pp
        Gk
        SSClAcc
        
        % Visualization Class
        Presenter
    end
    
    methods
        function obj = DataClassifier(p,c,f,pres)
            %DATACLASSIFIER Construct an instance of this class (l,p)
            %instances of dataloader and dataprocessor
            obj.processor       = p;
            %Parameter for classifier training
            obj.SelChans        = c;
            obj.SelFreqs        = f;
            obj.NumSelFeatures  = length(c);
            obj.Ck              = [];
            obj.NumWins         = [];
            obj.F               = [];
            obj.U               = [];
            obj.P               = [];
            obj.fullfreqs       = [];
            obj.LabelIdx        = [];
            
            obj.CkTotal         = [];
            obj.NumWinsTotal    = [];
            obj.FTotal          = [];
            obj.UTotal          = [];
            obj.PTotal          = [];
            obj.fullfreqsTotal  = [];
            
            % Accuracy member initialization
            obj.SSAcc               = []; 
            obj.pp                  = [];
            obj.Gk                  = [];
            obj.SSClAcc             = [];
            
            obj.Presenter       = pres;
            [~, obj.SelChansId] = ismember(obj.SelChans, obj.processor.loader.channelLb);
            
            obj = obj.loadFromProcessor();
        end
        
        function obj = loadFromProcessor(obj)
            obj = obj.datasetCreator();
            obj = obj.createModel();
            obj = obj.computeTrainsetAccuracy();
            obj = obj.saveClassifier();
        end
        
        function obj = datasetCreator(obj)
            % Chain Ck from all session in DataProcessor
            for sId=1:obj.processor.loader.nsessions    
                
                
                if (obj.processor.loader.offlineRuns{sId} == 0)
                    fprintf("No offline to visualize for session %d\n",sId);
                else
                    obj.CkTotal = cat(1,obj.CkTotal,obj.processor.Ck{sId});
                    obj.NumWinsTotal = obj.NumWinsTotal + obj.processor.NumWins{sId};

                    obj.Ck{sId} = obj.processor.Ck{sId};
                    obj.NumWins{sId} = obj.processor.NumWins{sId};
                    
                    
                    obj.P{sId} = obj.processor.loader.sessionsDataOffline{sId}.P;
                    obj.fullfreqs{sId} = obj.processor.loader.sessionsDataOffline{sId}.freqs;
                    obj.PTotal = cat(1,obj.PTotal,obj.processor.loader.sessionsDataOffline{sId}.P);
                    obj.fullfreqsTotal = cat(1,obj.fullfreqsTotal,obj.processor.loader.sessionsDataOffline{sId}.freqs);
                    
                    luOLD  = length(obj.processor.U{sId});
                    luNEW = length(obj.fullfreqs);
                    
                    
                    %[~, idfreqs] = intersect(obj.fullfreqs, obj.processor.SelFreqs);
            
                    %obj.U = log(obj.P(:, idfreqs, :));
                    
                    
                    
                    SelFreqsId=obj.SelFreqs;
                    obj.FTotal = nan(obj.NumWinsTotal, obj.NumSelFeatures);
                    obj.F{sId} = nan(obj.NumWins{sId}, obj.NumSelFeatures);
                    for ftId = 1:obj.NumSelFeatures
                        cfrq  = SelFreqsId(ftId);
                        cchan = obj.SelChansId(ftId);
                        obj.F{sId}(:, ftId) = obj.processor.U{sId}(:, cfrq, cchan);
                    end

                    % Create dataset LabelIdx
                    obj.LabelIdx{sId} = obj.Ck{sId} == 771 | obj.Ck{sId} == 773;
                
                
                end
                
                
                
            end
            
            
        end
        
        function obj = createModel(obj)
            for sId=1:obj.processor.loader.nsessions 
                if (obj.processor.loader.offlineRuns{sId} == 0)
                    fprintf("No offline to visualize for session %d\n",sId);
                else
                    fprintf("Training model %d\n",sId);
                    obj.Models{sId} = fitcdiscr(obj.F{sId}(obj.LabelIdx{sId}, :), obj.Ck{sId}(obj.LabelIdx{sId}), 'DiscrimType','quadratic');
                    obj.Presenter.PresentClassifier(obj.F{sId},obj.Ck{sId},obj.Models{sId},obj.LabelIdx{sId},obj.SelChans,obj.SelFreqs)
                end
            end
            %obj.Model = fitcdiscr(obj.F(obj.LabelIdx, :), obj.Ck(obj.LabelIdx), 'DiscrimType','quadratic');
        end
        
        function obj = computeTrainsetAccuracy(obj)
            for sId=1:obj.processor.loader.nsessions 
                if (obj.processor.loader.offlineRuns{sId} == 0)
                    fprintf("No offline to visualize for session %d\n",sId);
                else
                    fprintf("Computing accuracy model %d\n",sId);
                    NumClasses = length(obj.processor.loader.classId);
                    [obj.Gk{sId}, obj.pp{sId}] = predict(obj.Models{sId}, obj.F{sId});

                    obj.SSAcc{sId} = 100*sum(obj.Gk{sId}(obj.LabelIdx{sId}) == obj.Ck{sId}(obj.LabelIdx{sId}))./length(obj.Gk{sId}(obj.LabelIdx{sId}));

                    obj.SSClAcc{sId} = nan(NumClasses, 1);
                    for cId = 1:NumClasses
                        cindex = obj.Ck{sId} == obj.processor.loader.classId(cId);
                        obj.SSClAcc{sId}(cId) = 100*sum(obj.Gk{sId}(cindex) == obj.Ck{sId}(cindex))./length(obj.Gk{sId}(cindex));
                    end
                end
                
            end
        end
        
        function obj = saveClassifier(obj)
            for sId=1:obj.processor.loader.nsessions 
                save(sId, 'obj.Model', 'obj.SelChansId', 'obj.SelFreqsId');
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

