%% DATA PRESENTER FUNCTION CLASS 
%{
Description: Collection of utils and functions userd for show 
              result of single data processing or single function

Authors:  Agostini Francesco (francesco.agostini.5@studenti.unipd.it)
          Deschaux Ophélie   (opheliecandicemarine.deschaux@studenti.unipd.it)
          Marcon Francesco   (francesco.marcon.2@studenti.unipd.it)

Version: 0.1

%}
%% Class DataPresenter definition
classdef DataPresenter
    %% Class Properties 
    properties
        version
        data
        % Display Settings
        display_input_data           
        display_erd_ers              
        display_fisher_score         
        display_classifier           
        display_accumulated_evidence 
    end
    
    %% Class constructor and loader method
    methods
        %% Constructor
        function obj = DataPresenter(i,e,f,c,a)
            obj.version                      = 1.0;
            obj.display_input_data           = i;
            obj.display_erd_ers              = e;
            obj.display_fisher_score         = f;
            obj.display_classifier           = c;
            obj.display_accumulated_evidence = a;
        end
        
        %% Function Present RawData
        % Static Function for raw data plot
        function obj = PresentRawData(obj,name,P,TYP,DUR,POS,Rk,Mk)
            if (obj.display_input_data == 0)
                return;
            end
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
            if (obj.display_erd_ers == 0)
                % nothing to display
            else
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
        end

        %% Fisher Score plot
        % for use comment loadSessions() in the constructor DataLoader()
        
        %     obj: is costructor object
        %     folder: is a array of select folder
        function obj = PresentFisherScore( obj, name, NumRuns,NumFreqs,NumChans,channelLb,freqs ,FisherScore)
            if (obj.display_fisher_score == 0)
                % nothing to display
            else
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
        end
        
        
        function obj = PresentClassifier(obj,name,F,Ck,Model,LabelIdx,SelChans,SelFreqs)
            if (obj.display_classifier == 0)
                % nothing to display
            else
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
                sgtitle(name);
            end
            
        end  
        
        
        function obj = PresentAccumulatedEvidence(obj,name,Tk,Ck,pp,ipp,NumTrials)
            if (obj.display_accumulated_evidence == 0)
                % nothing to display
            else
                %% Plot accumulated evidence and raw probabilities
                fig1 = figure;

                CueClasses    = [771 783 773];
                LbClasses     = {'both feet', 'rest', 'both hands'};
                ValueClasses  = [1 0.5 0];
                Threshold     = 0.7;

                %SelTrial = 50;
                SelTrial = 5;

                % Trial 15: good rest
                % Trial 80: bad rest
                % Trial 55: good both hands
                % Trial 58: good both feet

                cindex = Tk == SelTrial;

                [~, ClassIdx] = ismember(unique(Ck(cindex)), CueClasses);

                GreyColor = [150 150 150]/255;
                LineColors = {'b', 'g', 'r'};

                hold on;
                % Plotting raw probabilities
                plot(pp(cindex, 1), 'o', 'Color', GreyColor);

                % Plotting accumulutated evidence
                plot(ipp(cindex), 'k', 'LineWidth', 2);

                % Plotting actual target class
                yline(ValueClasses(ClassIdx), LineColors{ClassIdx}, 'LineWidth', 5);

                % Plotting 0.5 line
                yline(0.5, '--k');

                % Plotting thresholds
                yline(Threshold, 'k', 'Th_{1}');
                yline(1-Threshold, 'k', 'Th_{2}');
                hold off;

                grid on;
                ylim([0 1]);
                xlim([1 sum(cindex)]);
                legend('raw prob', 'integrated prob');
                ylabel('probability/control')
                xlabel('sample');
                title(['Trial ' num2str(SelTrial) '/' num2str(NumTrials) ' - Class ' LbClasses{ClassIdx} ' (' num2str(CueClasses(ClassIdx)) ')']);

                    %sgtitle(name);
            end
            
        end
        
        
        
    
    end % methods (Static)
end % DataLoader