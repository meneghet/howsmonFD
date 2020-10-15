function fig_handle = pROC(SE, FPday, parValues, varargin)

p = inputParser;
addParameter(p,'fig_h',[])
addParameter(p,'plotBest',0)
addParameter(p,'msize',8)
addParameter(p,'labels',{'','',''})
addParameter(p,'connect',0)
addParameter(p,'fontSize',17)
parse(p,varargin{:});
fig_handle = p.Results.fig_h;
plotBest = p.Results.plotBest;
msize = p.Results.msize;
connectMarkers = p.Results.connect;
labels = p.Results.labels;
fontSize = p.Results.fontSize;

%% SETTINGS
if isempty(fig_handle)
    fig_handle=figure('Color','w');
else
    fig_handle=figure(fig_handle);
end
% font size
set(gca,'FontSize',fontSize)

%% read data

% different colors for each unique value of the selected parameter
par1vals = parValues(:,1);
par1_unique = unique(par1vals);
my_cmap = parula(length(par1_unique));

%% PLOT

% Z is used to save a number for each dot to use custom marker callback
Z = 1:length(SE);

% for each value of chosen parameter
hold on
for k = 1:length(par1_unique)
    V1 = par1_unique(k);
    % select data with this parameter value
    selected = (par1vals == V1);
    % assign a color
    par_color = my_cmap(par1_unique == V1,:);
    % plot
    h(k) = plot3(SE(selected), FPday(selected), Z(selected),'o-');
    % set style
    set(h(k),'MarkerSize',msize,'MarkerFaceColor',par_color,'MarkerEdgeColor',par_color,'Color',par_color);
    if ~connectMarkers
        set(h(k),'Linestyle','none');
    end
end

%% if second parameter is specified
if connectMarkers
    if size(parValues,2) > 1
        % read second parameter values
        par2vals = parValues(:,2);
        par2_unique = unique(par2vals);
        % delete lines
        set(h,'Linestyle','none');
        
        for n1 = 1:length(unique(par1vals))
            V1 = par1_unique(n1);
            for n2 = 1:length(unique(par2vals))
                V2 = par2_unique(n2);
                
                % color based on parameter 1
                par_color = my_cmap(par1_unique == V1,:);
                % select data with same parameter 1 and parameter 2 value
                selected = (par1vals == V1) & (par2vals == V2);
                % connect markers on this data
                plot3(SE(selected), FPday(selected), Z(selected),'-','Color',par_color);
                
            end
        end
    end
end



%% axis labels
if isempty(labels{1})
    xlabel('SE')
else
    xlabel(labels{1});
end
if isempty(labels{2})
    ylabel('FP/day')
else
    ylabel(labels{2});
end

axis([0 1 0 1]);
grid on
box on

%% colorbar (function below)
cb = make_my_colorbar(my_cmap, par1_unique, labels{3}, fontSize);

%% best setting
if plotBest
    
    distance = (1-SE).^2 + FPday.^2;
    [a,b] = min(distance);
    
    x = SE(b);
    y = FPday(b);
    
    h = plot(x,y,'o','MarkerSize',15,'Color','k','LineWidth',2);
    legend(h,['best'],'Location','northwest');
end

end

function cb = make_my_colorbar(my_cmap, par_list, label, fontsize)

colormap(my_cmap);
cb = colorbar;
caxis([min(par_list) max(par_list)])

if length(par_list)>12
    step = 2;
    set(cb,'Ticks',par_list(1:step:end))
else
    step = (par_list(end)-par_list(1))/length(par_list);
    set(cb,'Ticks',step/2+par_list(1):step:par_list(end)-step/2)
    set(cb,'TickLabels',arrayfun(@num2str, par_list, 'UniformOutput', false))
end

cb_title = get(cb,'Title');
set(cb_title,'String',label,'FontSize',fontsize,'FontWeight','bold')

end