function h = plotheatmap(groupNames, p_matrix)
    
    p_matrix(triu(true(size(p_matrix)), 1)) = NaN;
    ps4plotting =  p_matrix(2:end,1:end-1);
    figure('Theme','light'); 
    h = heatmap(groupNames(1:end-1), groupNames(2:end), ps4plotting, ...
        'ColorLimits', [0 0.05], ...
        'Colormap', flipud(hot), ...
        'CellLabelFormat','%.2e');

    h.MissingDataLabel = '';
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    

    % Annotate significance
    mask = ps4plotting < 0.05;
    p = ps4plotting(mask);
    [r, c] = find(ps4plotting < 0.05);
    ax = struct(h).Axes;  % access underlying axes of heatmap
    for k = 1:length(r)
        x = c(k);
        y = r(k);
        stars = siglevel(p(k));
        % Get center of each cell in axis units
        text(ax, x-.1, y-.5, stars, 'HorizontalAlignment', 'right', ...
             'VerticalAlignment', 'top', 'FontSize', 30, 'Color', 'b', 'FontWeight','bold');
    end
    ax = gca; ax.FontSize = 20; 

end

function stars = siglevel(p)

    if p<=1e-3
        stars='***'; 
    elseif p<=1e-2
        stars='**';
    elseif p<=0.05
        stars='*';
    else
        stars='';
    end
end