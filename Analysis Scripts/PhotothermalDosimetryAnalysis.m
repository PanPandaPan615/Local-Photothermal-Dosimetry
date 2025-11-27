%% Load Data
clc, clear, close all
cd ..,  addpath(genpath(cd))

%Single Pulse Diode
PW = load('TheoreticalBeamWidthData.mat');
PW = PW.PW;
Control = PW(:,end);

% Normalize Traces
for i = 1:length(PW)
    PW(i).Normalize = RemoveDC(PW(i).MeanF,1500, 250, 1);
end

% Constants
t = ((0:length(Control.MeanF)-1)*0.02-30)';
GreenCM = [zeros(256, 1),linspace(0, .75, 256)', zeros(256, 1)];
Titles = {'350 \mus','2 ms', '8 ms','Control'};
res = 0.288798; % Optical Resolution um/px
%%
halfwidth = 241; % measured halfwidth
% halfwidth = SpotSizeCalc(0.22,1,1.51,400,200,100); % theoretical halfwidth
SpotSize = (pi*halfwidth^2)*1e-8;
Fs = 50;

lambda = 1470;
mua = linterp(28.4,21.23,1460,1480,lambda);%  interpolation of Hale and Querry 1974 (cm^-1)
density_water = .9932; % density of water at 37C
cp = 4.18; % isobaric specific heat kJ/(kg*K) from engineering toolbox

LOOP = 1:length(PW);
%H50comp = readmatrix('StimH50Comps.csv');

% Control 
preidx = 1500-250:1500;   postidx = 1501:(1501+1000);
Control.pre = mean(Control.Normalize(preidx,:),1);
Control.post = mean(Control.Normalize(postidx,:),1);
Control.dF = (Control.post - Control.pre);
Control.dFstats = datastats(Control.dF');
Control.bound = Control.dFstats.max;

%% Strength Response Curves
RErange = [8 10 12];     % Radiant Exposure Range to fit over
RateParameter = [];

figure(Theme='light'); set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
for i = LOOP(1:3)
    binstep = iqr(PW(i).RE)*2/nthroot(length(PW(i).RE),3); % Freedman–Diaconis rule

    pre = mean(PW(i).Normalize(preidx,:),1);
    post = mean(PW(i).Normalize(postidx,:),1);
    PW(i).dF = (post-pre);
    
    RE = 0.1:0.1:RErange(i);
    PWColors(i,:) = [i*0.1/2 i*0.25 i*0.1/2];
    PW(i).RErange = RE;
    
    [PW(i).fit, PW(i).gof, PW(i).MeandF, PW(i).dFerr,PW(i).MeanRE] =...
        MeanResponse(PW(i),1500,1000,binstep);

    PW(i).ResponseCurve = feval(PW(i).fit,RE);
    PW(i).CI = predint(PW(i).fit, RE, 0.95,'functional','off');
    plot(RE,PW(i).ResponseCurve,LineWidth=2,Color=PWColors(i,:)); hold on
    errorbar(PW(i).MeanRE, PW(i).MeandF, PW(i).dFerr, 'vertical', 'LineStyle', 'none', 'Color', 'k', 'CapSize', 5,'HandleVisibility','off');
    scatter(PW(i).MeanRE,PW(i).MeandF,50,'filled','MarkerEdgeColor','k','MarkerFaceColor',PWColors(i,:),'HandleVisibility','off')
   
    conf = confint(PW(i).fit, 0.95); % 95% confidence intervals
    SE = (conf(2, 2)-conf(1, 2))/2; % standard error of b
    SE = [0 ,SE, -SE];
    RateParameter(i,:) = PW(i).fit.b + SE;

    PW(i).gof.rsquare = round(PW(i).gof.rsquare,3);
    PW(i).AUC = sum(PW(i).Normalize(postidx,:));
end

xlim([0 12]), ylim([-.1 7]), hold off
xlabel('Radiant Exposure (J/cm^2)');   ylabel('Mean \DeltaF/F (A.U.)');
set(gca, 'FontSize',24);
set(gca, 'XTick', 0:2:20, 'XTickLabels', 0:2:20);
legend('350 \mus', '2 ms', '8 ms',...
        "Location", "NorthWest","Color",'none','EdgeColor','none');

labels = repmat(1:3, [3 1])';
figure(Theme='light'); 
    b = bar(RateParameter(:,1)); hold on; % excluding pulse train because the fit was so bad
    b.FaceColor = 'flat';         
    for i = 1:numel(RateParameter(:,1)), b.CData(i,:) = PWColors(i,:); end
    errorbar(1:3, RateParameter(:,1), RateParameter(:,2)-RateParameter(:,1), 'k.', 'LineWidth', 1.5);
    ylabel('Exponential Growth Rate');
    set(gca, 'XTick', 1:3, 'XTickLabels', Titles(1:3));
    set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
    set(gca,'Fontsize', 24)

    [~,~,stats] = anova1(RateParameter(:),labels(:),"off");
    c = multcompare(stats,'Display','off','CriticalValueType','bonferroni');
    p =  c(:, 6); mask = p<0.05; 
    pairs =  num2cell(c(mask, 1:2),2);
    sigstar(pairs, p(mask)); hold off

%% % Grouping & Mean Reaction for each pulse width
for i = 1:3
    PW(i).MaskUnreactive = PW(i).dF<Control.bound; % Unreactive Mask
    PW(i).MaskHigh = PW(i).dF>1; % High Mask
    PW(i).MaskLow = ~PW(i).MaskUnreactive & ~PW(i).MaskHigh; % Middle Mask

    temp = vertcat(PW(i).MaskUnreactive , PW(i).MaskLow , PW(i).MaskHigh );

    figure(Theme='light');   set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 5]);
    ClassPlot(PW(i), temp, PWColors, Titles(i),Fs);  
    set(gca, 'FontSize',24);  xlim([0 90])

    PW(i).MaskReactive = PW(i).MaskLow|PW(i).MaskHigh;
    PW(i).dFStats = datastats(PW(i).dF(PW(i).MaskReactive)');
    PW(i).REStats = datastats(PW(i).RE(PW(i).MaskReactive)');
end

% Plotting individual traces in each class
% close all
% for i = LOOP(1:3)
%    t = (1:length(PW(i).Normalize(:,PW(i).MaskUnreactive)))*0.02-30;
%    figure; 
%    subplot(3,1,1),   plot(t,PW(i).Normalize(:,PW(i).MaskUnreactive));
%    title('Unreactive'), ylabel('\DeltaF/F (A.U.)');
%    xlim([-30 60]), ax = gca; ax.FontSize = 24;
%    set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 5]);
% 
%    subplot(3,1,2),   plot(t,PW(i).Normalize(:,PW(i).MaskLow));
%    title('Low Amplitude'), ylabel('\DeltaF/F (A.U.)');
%    xlim([-30 60]), ax = gca; ax.FontSize = 24;
%     set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 5]);
% 
%    subplot(3,1,3),   plot(t,PW(i).Normalize(:,PW(i).MaskHigh)');
%    title('High Amplitude'), ylabel('\DeltaF/F (A.U.)'), xlabel("Time (sec)")
%    xlim([-30 60]), ax = gca; ax.FontSize = 24;
%    set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 10]);
% 
%    temp = strjoin(['Pulse Width = ' Titles(i)]);
%    sgtitle(temp,'FontSize',36)
% end

%% Probit Plots
RErange = 0:0.01:13;        ActivationProb = 0:0.01:.9;

figure(Theme='light');
set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
for i = LOOP(1:3)
    binsize = iqr(PW(i).RE)*2/nthroot(length(PW(i).RE),3); % Freedman-Diaconis
    mask = PW(i).MaskReactive; %sum(mask)
    
    % Probit Creation
    [PW(i).Probit, PW(i).bins, PW(i).cdf, PW(i).Probitgof] =...
        ProbitAnalysis(PW(i).RE,mask,binsize,0);
    PW(i).ProbitCurve = feval(PW(i).Probit,RErange);
    PW(i).ProbitCI = predint(PW(i).Probit, RErange,0.95,'functional');
    Probit = [PW(i).ProbitCurve PW(i).ProbitCI];
    conf = confint(PW(i).Probit, 0.95); % 95% confidence intervals 
    SE = (conf(2, 1)-conf(1, 1))/2; % standard error of H50
    SE = [0 ,SE, -SE];
    H50(i,:) = PW(i).Probit.a+SE;

    mpeak(i,:)  = max(diff(Probit,1,1),[],1)/0.01; 
    
    % Group Plotting
    plot(RErange,PW(i).ProbitCurve,LineWidth=3,Color=PWColors(i,:)); hold on 
    scatter(PW(i).bins,PW(i).cdf,50,'filled','MarkerEdgeColor','k','MarkerFaceColor',PWColors(i,:),'HandleVisibility','off')
   
end
xlabel('Radiant Exposure (J/cm^2)'), ylabel('Activation Probability')
legend(Titles(1:3),"Location","northwest","Color",'none','EdgeColor','none')
ylim([0 1.1]), xlim([0 12]),
set(gca,'Fontsize', 24)

labels = repmat(1:3, [3 1])';

% Probit Slopes
figure(Theme='light');
    set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
    b = bar(mpeak(:,1));  hold on;
    b.FaceColor = 'flat';
    for i = 1:numel(mpeak(:,1)), b.CData(i,:) = PWColors(i,:); end 
    errorbar(1:numel(mpeak(:,1)), mpeak(:,1), mpeak(:,2)-mpeak(:,1),  mpeak(:,3)-mpeak(:,1), 'k.', 'LineWidth', 1.5);
    title("Slope of Probit Curves"), ylabel('m_p_e_a_k (J/cm^2)^-^1');
    set(gca, 'XTick', 1:length(Titles(1:3)), 'XTickLabels', Titles(1:3));
    set(gca, 'FontSize',24);

    [~,~,stats] = anova1(mpeak(:),labels(:),"off");
    c = multcompare(stats,'Display','off','CriticalValueType','bonferroni');
    p =  c(:, 6); mask = p<0.05;
    pairs =  num2cell(c(mask, 1:2),2);
    sigstar(pairs, p(mask)); hold off
    
   
% H50
figure(Theme='light');
    set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
    b = bar(H50(:,1)'); hold on;
    b.FaceColor = 'flat';
    for i = 1:numel(mpeak(:,1)), b.CData(i,:) = PWColors(i,:); end 
    errorbar(1:numel(H50(:,1)), H50(:,1), H50(:,2)-H50(:,1),  H50(:,3)-H50(:,1), 'k.', 'LineWidth', 1.5);
    set(gca, 'XTick', 1:length(Titles(1:3)), 'XTickLabels', Titles(1:3));
    ylabel('H_5_0 (J/cm^2)'); ylim([0 15])
    set(gca,'Fontsize', 24)

    [~,~,stats] = anova1(H50(:),labels(:),"off");
    c = multcompare(stats,'Display','off','CriticalValueType','bonferroni');
    p =  c(:, 6); mask = p<0.05;
    pairs =  num2cell(c(mask, 1:2),2);
    sigstar(pairs, p(mask)); hold off

%% Polar Plots
meddist = []; maddist = [];  reactingdistances = [];
RErange = 0:0.01:13;  Threshlevel = 0.01;
for i = LOOP(1:3)
    [~, idx] = min(abs(PW(i).ProbitCurve-Threshlevel));
    mask = PW(i).MaskReactive & (PW(i).RE>RErange(idx));
    
    figure(Theme='light');
    set(gcf, 'Units', 'inches', 'Position', [1, 1, 10, 5]);
        polarscatter(PW(i).Theta, PW(i).Rho*res, 20, PW(i).RE,'filled'),  colormap("turbo");
        cb = colorbar("AxisLocation","out","FontSize",24,"Location","eastoutside"); clim([0 10])
        title(Titles(i)), ylabel(cb, 'J/cm^2',"FontSize",24); 
        set(gca,'FontSize',24); ax.ThetaTickLabel = {}; rlim([0 500]);
        grid on
    
        meddist(i) = median(PW(i).Rho(mask))*res;
        maddist(i)  = mad(PW(i).Rho(mask))*res;
        reactingdistances = [reactingdistances PW(i).Rho(mask)*res ];
        RadialThresh(i) =  meddist(i)+maddist(i)*3;
        Rho_RE_corr(i) = corr(PW(i).Rho'.*res,PW(i).RE',"Type","Spearman");
end

%% Reacting distributions
ticks = concatCellVectors({'U','L', 'H'},{Titles{1:3}});
figure(Theme='light');  set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 5]);
RE = []; glabels = []; pwlabels = []; pairs = []; p =[];
k = 0;
  for i = LOOP(1:3)
    mask = PW(i).MaskUnreactive & PW(i).Rho*res<RadialThresh(i);
    yU = PW(i).RE(mask);  xU = ones(1,length(yU))+k;
    mask = PW(i).MaskLow & PW(i).Rho*res<RadialThresh(i);
    yL = PW(i).RE(mask);         xL = ones(1,length(yL))*2+k;
    mask = PW(i).MaskHigh & PW(i).Rho*res<RadialThresh(i);
    yH = PW(i).RE(mask);        xH = ones(1, length(yH))*3+k;

    swarmchart(xU, yU, 20,"filled",'MarkerEdgeColor',...
              'k',MarkerFaceColor=PWColors(1,:),XJitterWidth=0.7); hold on
    swarmchart(xL, yL, 20,"filled",'MarkerEdgeColor',...
              'k',MarkerFaceColor=PWColors(2,:),XJitterWidth=0.3);
    swarmchart(xH, yH, 20,"filled",'MarkerEdgeColor',...
              'k',MarkerFaceColor=PWColors(3,:),XJitterWidth=0.3);
        
    [~,tbl,stats] = kruskalwallis( [yU, yL, yH] ,[xU, xL, xH],'off');
    c = multcompare(stats,'Display','off',CriticalValueType='dunn-sidak');
    pairs = [pairs; num2cell(c(:, 1:2)+k,2)];
    p = [p; c(:, 6)];

    k = k+3;
  end

ylabel('Radiant Exposure (J/cm^2)'); ylim([0 11])
set(gca, 'XTick', [2,5,8], 'XTickLabels', Titles(1:3), 'FontSize',24);
mask = p<0.05;
sigstar(pairs(mask), p(mask)); hold off;
legend({'Unreactive', 'Low Amplitude', 'High Amplitude'},...
        "Location", "southeast", 'Orientation', 'horizontal',...
        "Color",'none','EdgeColor','none'); % Assuming Titles has labels for each subplot

p = reshape(p,[i i]);
figure(Theme='light')
heatmap({Titles{1:3}},{'U-L','U-H','L-H'},p,'ColorLimits', [0 0.05], ...
        'Colormap', flipud(hot), ...
        'CellLabelFormat','%.2e')
        set(gca, 'FontSize',20);
%% Mean Radiant Exposure of Reacting Cells
Threshlevel = 0.01;
LocalRE = []; SpotRE = []; groupLabels=[]; REmeans = []; REdevs = []; REerr=[]; DeltaTemp = [];
    for i = LOOP(1:3)
        mask = PW(i).MaskReactive;         sum(mask);
        [~, idx] = min(abs(PW(i).ProbitCurve-Threshlevel));
        PW(i).LocalRE =   PW(i).RE(mask & (PW(i).RE>RErange(idx)));
        PW(i).SpotRE =  PW(i).PulseEnergy(mask & (PW(i).RE>RErange(idx)))/SpotSize*1e-3;
                        
        LocalRE = [LocalRE  PW(i).LocalRE]; SpotRE = [SpotRE  PW(i).SpotRE];
        groupLabels = [groupLabels; ones(length(PW(i).LocalRE),1)*i];

        temp1 = [mean(PW(i).SpotRE); mean(PW(i).LocalRE); ];   REmeans = [REmeans temp1];
        temp = [std(PW(i).SpotRE); std(PW(i).LocalRE); ];     REdevs = [REdevs temp];
        temp = temp./sqrt(length(PW(i).LocalRE));             REerr = [REerr temp];

        PW(i).Temp = (mua*PW(i).LocalRE)./(density_water*cp); % Temperature conversion;
    end

groupLabels = [groupLabels; groupLabels+length(LOOP(1:3))];
data = [SpotRE LocalRE]';

% Reacting distributions
figure(Theme='light');  set(gcf, 'Units', 'inches', 'Position', [2, 2, 14, 7]);
temp = [1 2 3 1 2 3];
  for i = unique(groupLabels)'
      mask = groupLabels == i;
      sum(mask);
      x = ones(1, sum(mask))*i;
      swarmchart(x, data(mask), 50,"filled",'MarkerEdgeColor',...
       'k',MarkerFaceColor=PWColors(temp(i),:),XJitterWidth=0.3); hold on
  end
    
  % Intragroup comps
    for i = 1:2
        mask = logical(sum(groupLabels==((1:3)+(i-1)*3),2));
        [~,~,statsVal] = anova1(data(mask), groupLabels(mask),'off'); 
        c = multcompare(statsVal,'CriticalValueType','bonferroni','Display','off');
        c(:,1:2) = c(:,1:2)+(i-1)*3;
        pmask = c(:,6)<0.05;
        pairs = num2cell(c(pmask, 1:2),2);
        sigstar(pairs, c(pmask,6)); 
    end

    % Intergroup comps
    for i = 1:3
        [~,p] = ttest(data(groupLabels==i), data(groupLabels==(i+3)),"Tail","both");
        pairs = [i, i+3];
        pairs = num2cell(pairs,2);
        sigstar(pairs, p); 
    end
   
    ylabel('Radiant Exposure (J/cm^2)');
    set(gca, 'XTick', [2,5], 'XTickLabels', {'Average Dosage' 'Local Dosage'});
    set(gca, 'FontSize',24);
    legend(Titles(1:3),"Location", "southeast", "Color",'none',...
        'EdgeColor','none','Orientation','horizontal'); % Assuming Titles has labels for each subplot
    ylim([0 14]), xlim([0.5 6.5]),    hold off;
    
%% Rad Exposure Comparisons
meanH50 = reshape(H50comp(:,1),[3,2])';
H50err = reshape(H50comp(:,2),[3,2])'-meanH50;
figure(Theme='light'),   set(gcf, 'Units', 'inches', 'Position', [2, 2, 10, 7]);
    b = bar(meanH50); colororder(PWColors), hold on
    errorbar([b.XEndPoints], meanH50(:), H50err(:), 'k.', 'LineWidth', 1.5);
    set(gca, 'XTickLabels', {'Theoretical','TLA'});
    ylabel('H_5_0 (J/cm^2)'); ylim([0 15])
    ax = gca; ax.FontSize = 24;
    groupLabels = repmat(sort([b.XEndPoints]), [numel(H50comp(1,:)) 1])';
    
    % Intragroup comps
    for i = 1:2
        mask = round(groupLabels)==i;
        [~,~,statsVal] = anova1(H50comp(mask), groupLabels(mask),'off'); 
        c = multcompare(statsVal,'Display','off');
        c(:,1:2) = c(:,1:2)+(i-1)*3;
        pmask = c(:,6)<0.05;
        pairs = [groupLabels(c(pmask,1),1), groupLabels(c(pmask,2),1)];
        pairs = num2cell(pairs(pmask,:),2);
        sigstar(pairs, c(pmask,6)); 
    end

    % Intergroup comps
    for i = 1:3
        [~,p] = ttest(H50comp(i,:), H50comp(i+3,:),"Tail","right");
        pairs = [groupLabels(i,1), groupLabels(i+3,1)];
        pairs = num2cell(pairs,2);
        sigstar(pairs, p); 
    end

    legend(Titles(1:3),"Location", "north", "Color",'none','EdgeColor','none',...
        'Orientation','horizontal'); 

%% Damage Calcium traces
DamageE = [7.2 7.5];
for i = 2:3
    mask = PW(i).PulseEnergy<DamageE(i-1) & PW(i).MaskReactive;
    t = (1:length(PW(i).Normalize(:,1)))*0.02-30;
    figure(Theme='light'),  set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 7]);
    plot(t,PW(i).Normalize(:,mask));
    ylabel('\DeltaF/F (A.U.)'); xlabel('Time (sec)');
    ylim([-0.1 4.5]), xlim([t(1) t(end)]),   grid on; 
    temp = strjoin(['Pulse Width = ' Titles(i)...
           '  Pulse Energy \leq ' num2str(DamageE(i-1)) 'mJ']);
    title(temp,'FontSize',36),  set(gca, 'FontSize',24);
end

%% FUNCTIONS
function interped = linterp(x1,x2,y1,y2,desired)
    m = (y2 - y1)/(x2-x1);
    interped = (desired-y1)/m+x1;
end


function  Data = AppendData(Data1,Data2,idx,Plot)
   

    Data.PulseEnergy =     cat(2,Data1.PulseEnergy,    Data2.PulseEnergy);
    Data.RE =     cat(2,Data1.RE,    Data2.RE);
    Data.Energy = cat(2,Data1.Energy,Data2.Energy);
    Data.Area =   cat(2,Data1.Area,  Data2.Area);

    Data.MeanF =   cat(2,Data1.MeanF(idx,:),  Data2.MeanF(idx,:));
    Data.Std =     cat(2,Data1.Std(idx,:),    Data2.Std(idx,:));
    Data.MaxF =    cat(2,Data1.MaxF(idx,:),   Data2.MaxF(idx,:));
    Data.Entropy = cat(2,Data1.Entropy(idx,:),Data2.Entropy(idx,:));
    Data.Skw =     cat(2,Data1.Skw(idx,:),    Data2.Skw(idx,:));
    Data.Uni =     cat(2,Data1.Uni(idx,:),    Data2.Uni(idx,:));

    Data.Centroids =    cat(2,Data1.Centroids,  Data2.Centroids);
    Data.Rho =          cat(2,Data1.Rho,  Data2.Rho);
    Data.Theta =        cat(2,Data1.Theta,  Data2.Theta);

    if Plot == 1
        figure;
        subplot(2,3,1);
        temp =  RemoveDC(Data.MeanF,500,1);
        plot(temp)
        title('Raw Mean ROI Fluorescence'), xticklabels(-30:20:60);
        xlabel('Time (s)'), ylabel('Mean Gray Value \DeltaF/F (A.U.))')

        subplot(2,3,2);
        temp =  RemoveDC(Data.Std,500,1);
        plot(temp)
        title('ROI Fluorescence Standard Deviation'), xticklabels(-30:20:60);
        xlabel('Time (s)'), ylabel('\sigma Gray Value \DeltaF/F (A.U.)')

        subplot(2,3,3);
        temp =  RemoveDC(Data.MaxF,500,1);
        plot(temp)
        title('Raw Max ROI Fluorescence'),xticklabels(-30:20:60);
        xlabel('Time (s)'), ylabel('Max Gray Value \DeltaF/F (A.U.)')

        subplot(2,3,4);
        temp =  RemoveDC(Data.Entropy,500,0);
        plot(temp)
        title('ROI Entropy'), xticklabels(-30:20:60);
        xlabel('Time (s)'), ylabel('(A.U.)')

        subplot(2,3,5);
        temp =  RemoveDC(Data.Skw,500,0);
        plot(temp)
        title('ROI Skewness'), xticklabels(-30:20:60);
        xlabel('Time (s)'), ylabel('(A.U.)')

        subplot(2,3,6);
        temp =  RemoveDC(Data.Uni,500,0);
        plot(temp)
        title('ROI Uniformity'), xticklabels(-30:20:60);
        xlabel('Time (s)'), ylabel('(A.U.)')
    end

end

function  ClassPlot(PW, masks, PWColors, Title, Fs)
    for i = 1:size(masks,1)
        temp  = masks(i,:);
        MeanPlot(PW.Normalize(:,temp),Title,PWColors(i,:),Fs); hold on
    end

    legend('Unreactive','Low Amplitude','High Amplitude',...
           "Location", "northwest","Color",'none','EdgeColor','none');
end

function MeanPlot(traces,TITLE,COLOR,Fs)
    TraceMean = mean(traces,2);
    TraceSTD = std(traces,0,2);
    TraceErr = TraceSTD./sqrt(min(size(traces)));
    UB = TraceMean+TraceErr; LB=TraceMean-TraceErr;
    t = (1:length(TraceMean))./Fs;

    plot(t,TraceMean,LineWidth=1,Color=COLOR), hold on
    fill([t, fliplr(t)], [UB', fliplr(LB')], COLOR, 'FaceAlpha', 0.2, 'EdgeAlpha', 0,'HandleVisibility','off'); hold off
    title(TITLE), xlabel('Time (s)'), ylabel('\DeltaF/F (A.U.)'),     ax = gca; ax.FontSize = 20;
end

function ShadedErr(X, CI, Color, Alpha)

    fill([X, fliplr(X)], [CI(:, 1)', fliplr(CI(:, 2)')],...
     Color, 'FaceAlpha', Alpha, 'EdgeAlpha', 0,'HandleVisibility','off');

end

function [fitResult, gof, Delta, err,REVector] = MeanResponse(PW,centeridx,range,binstep)
    
    REVector = 0:binstep:max(PW.RE);
    idx = 1;
    for i = REVector
        mask = (i<PW.RE)&(PW.RE<i+binstep);
        n(idx) = sum(mask);
        pre = PW.MeanF(centeridx-range:centeridx-1,mask);    post = PW.MeanF(centeridx+1:end,mask);
        means = [mean(pre,"all") mean(post,"all")];
        stds = [std(pre,0,"all") std(post,0,"all")];
        err(idx) = sqrt(stds(1)^2+stds(2)^2)/n(idx)/means(1);
        Delta(idx) = (means(2)-means(1))/means(1);
        idx = idx+1;
    end
    
    % Define the weights as the inverse of the standard deviations
    weights = 1 ./ err;
    weights(isnan(weights))=0;
    %mask = Delta > fitthresh; for piecewise fit
    
    % Define the exponential model
    expModel = fittype('a*exp(b*x)+c');
    
    % Create a fitoptions object specifying the weights
    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                      'Weights', weights, ...
                      'StartPoint', [0, 1, 0]); % Initial guess for a and b
    
    % Fit the model to the data
    [fitResult, gof] = fit(REVector', Delta',expModel, opts);

end

function result = concatCellVectors(vec1, vec2)
    % Convert all elements to strings (if they are numeric)
    strVec1 = cellfun(@string, vec1, 'UniformOutput', false);
    strVec2 = cellfun(@string, vec2, 'UniformOutput', false);

    % Generate all combinations using ndgrid
    [A, B] = ndgrid(1:numel(strVec1), 1:numel(strVec2));

    % Concatenate strings with a space between
    result = arrayfun(@(i, j) strcat(strVec1{i}, " ", strVec2{j}), A(:), B(:), 'UniformOutput', false);
end

function [ShiftedTraces, trough] = RemoveDC(Traces,Stimpoint,MinFrame,Normalize)
ShiftedTraces = zeros(size(Traces));

    for i=1:min(size(Traces)) 
        trace = Traces(1:Stimpoint,i);
        smoothed = smooth(trace,5);
        trough = min(smoothed(end-MinFrame:end));
    
        ShiftedTraces(:,i) = Traces(:,i)-trough;
    
        if Normalize == 1
            ShiftedTraces(:,i) = ShiftedTraces(:,i)/trough;
        end
    end
end