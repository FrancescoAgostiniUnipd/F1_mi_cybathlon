%%  plot_Classifier(Model, F, LabelIdx, Ck, Channels, Frequences, choose)
%   Model,  Model
%   F,      X data insert into the model,
%   LabelIdx, value of classifier
%   Ck,      group of value
%   Channels, Channels selected
%   Frequences, Frequences selected,
%   choose,  dimension of the array to plot
function plot_Classifier(Model, F, LabelIdx, Ck, Channels, Frequences, choose)
    % fig2 = figure;
     gscatter(F(LabelIdx, choose(1)),F(LabelIdx, choose(2)),Ck(LabelIdx),'mg','.',[],'off');
    grid on;
    xlim([-8 0]);
    ylim([-8 1.5]);
    % xlim auto
    % ylim auto
    xlabel([Channels{choose(1)} '@' num2str(Frequences(choose(1))) 'Hz']);
    ylabel([Channels{choose(2)} '@' num2str(Frequences(choose(2))) 'Hz']);
    title(['Plot of ' Channels{choose(1)} '@' num2str(Frequences(choose(1))) 'Hz | ' Channels{choose(2)} '@' num2str(Frequences(choose(2))) 'Hz']);
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
    f = @(x1,x2) K + L(choose(1))*x1 + L(choose(2))*x2 + Q(choose(1),choose(1))*x1.^2 + ...
        (Q(choose(1),choose(2))+Q(choose(2),choose(1)))*x1.*x2 + Q(choose(2),choose(2))*x2.^2;
    
    h2 = fimplicit(f);
    h2.Color = 'r';
    h2.LineWidth = 2;
    h2.DisplayName = 'Boundary between boht hands & both feet';
    legend('both feet', 'both hands', 'Boundary');
    hold off;
end