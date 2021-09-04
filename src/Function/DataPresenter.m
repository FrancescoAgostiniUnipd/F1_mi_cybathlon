%% DATA PRESENTER FUNCTION CLASS 
%{
Description: Collection of utils and functions userd for show 
              result of single data processing or single function

Authors:  Agostini Francesco (francesco.agostini.5@studenti.unipd.it)
          Deschaux Oph�lie   (opheliecandicemarine.deschaux@studenti.unipd.it)
          Marcon Francesco   (francesco.marcon.2@studenti.unipd.it)

Version: 0.1

%}
%% Class DataPresenter definition
classdef DataPresenter
    %% Class Properties 
    properties
        version;
        data
    end
    
    %% Class constructor and loader method
    methods
        %% Constructor
        function obj = DataPresenter()
            obj.version = 0.1;
        end
        
        %% Function Present RawData
        % Static Function for raw data plot
        function obj = PresentRawData(obj,name,P,TYP,DUR,POS,Rk,Mk)
            figure;
            subplot(3,2,1);
            %plot(P);             % Plot samples
            title('samples')
            subplot(3,2,2);
            plot(TYP);           % Plot session TYP vector 
            title('TYP')
            subplot(3,2,3);
            plot(DUR);           % Plot session DUR vector
            title('DUR')
            subplot(3,2,4);
            plot(POS);           % Plot session POS vector 
            title('POS')
            subplot(3,2,5);
            plot(Rk);            % Plot session Rk 
            title('Rk')
            subplot(3,2,6);
            plot(Mk);            % Plot session Mk 
            title('Mk') 
            sgtitle(name);
        end
        
        %% Function Present ERD/ERS plot
        % static function for session data display
        function obj = PresentErdErs(obj, name, MinTrialDur, wshift, ChannelSelected,nclasses,ERD,classId,freqs,channelLb,classLb,tCk)
            figure;
            fprintf("MinTrialDur = %f | wshift = %f \n",MinTrialDur,wshift);
            t = linspace(0, MinTrialDur*wshift, MinTrialDur);
            ChannelSelected = [7 9 11]; 

            chandles = [];
            for cId = 1:nclasses

                climits = nan(2, length(ChannelSelected));
                for chId = 1:length(ChannelSelected)
                    subplot(2, 3, (cId - 1)*length(ChannelSelected) + chId);
                    cdata = mean(ERD(:, :, ChannelSelected(chId), tCk == classId(cId)), 4);
                    
                    
                    imagesc(t, freqs, cdata');
                    set(gca,'YDir','normal');
                    climits(:, chId) = get(gca, 'CLim');
                    chandles = cat(1, chandles, gca);
                    colormap(hot);
                    colorbar;
                    title(['Channel ' channelLb{ChannelSelected(chId)} ' | ' classLb{cId}]);
                    xlabel('Time [s]');
                    ylabel('Frequency [Hz]');
                    line([1 1],get(gca,'YLim'),'Color',[0 0 0])
                end

            end
            set(chandles, 'CLim', [min(min(climits)) max(max(climits))]);
            sgtitle(name);
        end

        %% Fisher Score plot
        % for use comment loadSessions() in the constructor DataLoader()
        
        %     obj: is costructor object
        %     folder: is a array of select folder
        function obj = PresentFisherScore( obj, name, NumRuns,NumFreqs,NumChans,channelLb,freqs ,FisherScore)
            OfflineRuns = 1:NumRuns;
            climits = [];
            handles = nan(NumRuns, 1);
            fig1 = figure;
            SelChans={};
            SelFreqs=[];
            FisherScoretemp=FisherScore;
            for rId = 1:length(OfflineRuns)
                subplot(1, length(OfflineRuns), rId);
                imagesc(FisherScore(:, :, OfflineRuns(rId))');

                % To select the freq and chan with the highest fisher score
                A=FisherScoretemp(:,:,OfflineRuns(rId));
                val = 10;
                while(val>0.9)

                    [val,idx] = max(A(:));
                    if val>0.9
                        [row,col] = ind2sub(size(A),idx);
                        A(row,col)=0;
                        FisherScoretemp(row,col,:)=0;
                        SelChans=cat(2,SelChans,channelLb(col));
                        SelFreqs=cat(2,SelFreqs,freqs(row));
                    end
                end

                axis square;
                set(gca, 'XTick', 1:NumFreqs);
                set(gca, 'XTickLabel', freqs);
                set(gca, 'YTick', 1:NumChans);
                set(gca, 'YTickLabel', channelLb);
                xtickangle(-90);

                title(['Calibration run ' num2str(OfflineRuns(rId))]);

                climits = cat(2, climits, get(gca, 'CLim'));
                handles(OfflineRuns(rId)) = gca;
            end


            set(handles, 'CLim', [min(min(climits)) max(max(climits))]);
            %txt = ['FISHER SCORE of ',name];
            sgtitle(name);
        end
        
        
        function obj = PresentClassifier(obj,F,Ck,Model,LabelIdx,SelChans,SelFreqs)
           %% Visualize classifier
            fig2 = figure;
            h1 = gscatter(F(LabelIdx, 1),F(LabelIdx, 2),Ck(LabelIdx),'kb','ov^',[],'off');
            grid on;
            xlim([-8 0]);
            ylim([-8 1.5]);
            xlabel([SelChans{1} '@' num2str(SelFreqs(1)) 'Hz']);
            ylabel([SelChans{2} '@' num2str(SelFreqs(2)) 'Hz']);
            axis square;
            hold on

            % Linear
            % K = Model.Coeffs(1,2).Const;
            % L = Model.Coeffs(1,2).Linear;
            % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;

            % Quadratic
            K = Model.Coeffs(1,2).Const;
            L = Model.Coeffs(1,2).Linear;
            Q = Model.Coeffs(1,2).Quadratic;
            f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
                (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;

            h2 = fimplicit(f);
            h2.Color = 'r';
            h2.LineWidth = 2;
            h2.DisplayName = 'Boundary between boht hands & both feet';
            legend('both feet', 'both hands', 'Boundary');
            hold off; 
            
            
        end    
        
        
        
    
    end % methods (Static)
end % DataLoader