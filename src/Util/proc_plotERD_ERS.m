%% Visualization ERD/ERS
%   To use proc_plotERD_ERS(data, [7,8,9], t,tCk,4:2:48, ERD)

%       -data: Object of DataLoader
%       -space: linearspace()
%       -ChannelSelected: three selected channels
%       -tCk: type of event (Both Hand 773 or Both Feet 771)
%       -freqs: frequency of visualization on chart
%       -ERD: data to visualize

    function  proc_plotERD_ERS(data, ChannelSelected, space,tCk,freqs, ERD)
        figure;
    
        chandles = [];
        for cId = 1:data.nclasses
            
            climits = nan(2, length(ChannelSelected));
            for chId = 1:length(ChannelSelected)
                subplot(2, 3, (cId - 1)*length(ChannelSelected) + chId);
                cdata = mean(ERD(:, :, ChannelSelected(chId), tCk == data.classId(cId)), 4);
                imagesc(space, freqs, cdata');
                set(gca,'YDir','normal');
                climits(:, chId) = get(gca, 'CLim');
                chandles = cat(1, chandles, gca);
                colormap(hot);
                colorbar;
                title(['Channel ' data.channelLb{ChannelSelected(chId)} ' | ' data.classLb{cId}]);
                xlabel('Time [s]');
                ylabel('Frequency [Hz]');
                line([1 1],get(gca,'YLim'),'Color',[0 0 0])
            end
            
        end
        set(chandles, 'CLim', [min(min(climits)) max(max(climits))]);
    end
