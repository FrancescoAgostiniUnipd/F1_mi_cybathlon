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
        Ck
        LabelIdx
        SelChansId
        NumWins
        F
        P
        U
        % Model members
        Model
        
        % Accuracy Member
        SSAcc
        pp
        Gk
        SSClAcc
        
    end
    
    methods
        function obj = DataClassifier(p,c,f)
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
           
            [~, obj.SelChansId] = ismember(obj.SelChans, obj.processor.loader.channelLb);
            
            %obj = obj.loadFromProcessor();
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
                obj.Ck = cat(1,obj.Ck,obj.processor.Ck{sId});
                obj.NumWins = obj.NumWins + obj.processor.NumWins{sId};
                obj.P = cat(1,obj.P,obj.processor.loader.sessionsDataOffline{sId}.P);
            end
            
            obj.U = log(obj.P(:, obj.processor.SelFreqs, :));
            SelFreqsId=obj.SelFreqs;
            obj.F = nan(obj.NumWins, obj.NumSelFeatures);
            for ftId = 1:obj.NumSelFeatures
                cfrq  = SelFreqsId(ftId);
                cchan = obj.SelChansId(ftId);
                obj.F(:, ftId) = obj.U(:, cfrq, cchan);
            end
            
            % Create dataset LabelIdx
            obj.LabelIdx = obj.Ck == 771 | obj.Ck == 773;
        end
        
        function obj = createModel(obj)
            obj.Model = fitcdiscr(obj.F(obj.LabelIdx, :), obj.Ck(obj.LabelIdx), 'DiscrimType','quadratic');
        end
        
        function obj = computeTrainsetAccuracy(obj)
            NumClasses = length(obj.processor.loader.classId);
            [obj.Gk, obj.pp] = predict(obj.Model, obj.F);

            obj.SSAcc = 100*sum(obj.Gk(obj.LabelIdx) == obj.Ck(obj.LabelIdx))./length(obj.Gk(obj.LabelIdx));
            
            obj.SSClAcc = nan(NumClasses, 1);
            for cId = 1:NumClasses
                cindex = obj.Ck == obj.processor.loader.classId(cId);
                obj.SSClAcc(cId) = 100*sum(obj.Gk(cindex) == obj.Ck(cindex))./length(obj.Gk(cindex));
            end
        end
        
        function obj = saveClassifier(obj,filename)
            save(filename, 'obj.Model', 'obj.SelChansId', 'obj.SelFreqsId');
        end
        
        function obj = loadClassifier(obj,filename)
            MDL = load(filename);
            obj.SelChansId = MDL.SelChansId;
            obj.SelFreqsId = MDL.SelFreqsId;
            obj.Model = MDL.Model;
        end
    end
end

