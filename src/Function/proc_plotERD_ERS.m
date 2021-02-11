%% Visualization ERD/ERS
%   To use:
% proc_plotERD_ERS(cdata, datas.channelLb{ChannelSelected(chId)}, t, session1.freqs, datas.classLb{cId});

%       -cdata: Data to Plot
%       -nameChannel: Name of the channel to reppresent
%       -space: linearspace()
%       -freqs: frequency of visualization on chart
%       -classLb: Name of class to represent


    function  [chandles, climits] = proc_plotERD_ERS(cdata, nameChannel, space,freqs,classLb)
        %tCk, ERD)
        % figure;
    
        chandles = [];
        % for cId = 1:data.nclasses
            
        %     climits = nan(2, length(ChannelSelected));
        %     for chId = 1:length(ChannelSelected)
        %         subplot(2, 3, (cId - 1)*length(ChannelSelected) + chId);
        % cdata = mean(ERD(:, :, ChannelSelected, tCk == data.classId(cId)), 4);
        imagesc(space, freqs, cdata');
        set(gca,'YDir','normal');
        climits = get(gca, 'CLim');
        chandles = cat(1, chandles, gca);
        colormap(hot);
        colorbar;
        title(['Channel ' nameChannel ' | ' classLb]);
        xlabel('Time [s]');
        ylabel('Frequency [Hz]');
        line([1 1],get(gca,'YLim'),'Color',[0 0 0])
        %     end
            
        % end
        % set(chandles, 'CLim', [min(min(climits)) max(max(climits))]);
    end
