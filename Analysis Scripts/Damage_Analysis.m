clc, clear, close all
cd ..,  addpath(genpath(cd))

% Load data
H50comp = readmatrix('DamageH50Comparison.csv');
Diode = load('SPDamage.mat');   Diode = Diode.Diode;

categories = {'Control_(_-_)', '350 \mus', '2 ms',...
              '8 ms','2 ms Control_(_+_)', '8 ms Control_(_+_)'};
Colors = [0,    0,    0;
          0.05, 0.25, 0.05;
          0.1,  0.5,  0.1;
          0.15, 0.75, 0.15;
          0.75, 0.1,  0.3;
          0.5,  0.1,  0.3;];
 
Appendend = [Diode.Measurements];

% Generate Experimental Masks  [NegControl 350us 2msTest 8msTest 2msPos 8msPos]
temp = [0 .35 2 8];
for i = 1:length(temp)
    if i<3
        ExpMask(i,:)   = [Appendend.PW] == temp(i);
    else 
        ExpMask(i,:)   = [Appendend.PW] == temp(i) &...
                         [Appendend.PulseEnergy] < 10;

        ExpMask(i+2,:) = [Appendend.PW] == temp(i) &...
                         [Appendend.PulseEnergy] > 10;
    end
end
    ExpMask(2,3) = 0;% Throw out bad data point in 350 us group

% Create experimental structure
for i = 1:length(ExpMask(:,1))
    Exp(i).RE = vertcat(Appendend(ExpMask(i,:)).RadiantEnergy);
    Exp(i).PrePIintensity = vertcat(Appendend(ExpMask(i,:)).PrePIintensity);
    Exp(i).PostPIintensity = vertcat(Appendend(ExpMask(i,:)).PostPIintensity);

    Exp(i).PrePICount = [Appendend(ExpMask(i,:)).PrePICount];
    Exp(i).CalcienCount = [Appendend(ExpMask(i,:)).CalcienCount];
    Exp(i).PostPICount = [Appendend(ExpMask(i,:)).PostPICount];

    Exp(i).DoseLabels = i*ones(1,length(Exp(i).RE));
    Exp(i).CountLabels = i*ones(1,length(Exp(i).PrePICount));

end

%% Create Groups
DiodeCountDifference = [Exp.PostPICount] - [Exp.PrePICount];
DiodeLDRatio = DiodeCountDifference./[Exp.CalcienCount]*100;

CountLabels = [Exp.CountLabels];
LengthVector = grpstats(CountLabels, CountLabels, 'numel');

for i = 1:length(Exp)
    Max(i) = max(Exp(i).RE);
end

halfwidth = 240;
SpotSize = (pi*halfwidth^2)*1e-8;
SpotRE = 14.3/SpotSize*1e-3;

% Bar Plots
% Absolute Difference
Mean = grpstats(DiodeCountDifference, CountLabels, 'mean');
StdDev = grpstats(DiodeCountDifference, CountLabels, 'std');
Err = StdDev./sqrt(LengthVector);

figure(Theme='light'); %set(gcf, 'Units', 'inches', 'Position', [1, 1, 10, 10]);
    b = bar(Mean); hold on; % excluding pulse train because the fit was so bad
    b.FaceColor = 'flat';         
    for i = 1:numel(Mean), b.CData(i,:) = Colors(i,:); end
    errorbar(1:numel(Mean), Mean, Err, 'k.', 'LineWidth', 1.5);
    title('Absolute Diode Damage Assay'),   ylabel('\Delta Dead Cells');
    set(gca, 'XTick', 1:length(categories), 'XTickLabels', categories);
    set(gca,'Fontsize', 24)

    [~,~,statsVal] = anova1(DiodeCountDifference, CountLabels,'off'); 
    c = multcompare(statsVal,'CriticalValueType','bonferroni','Display','off');
    p =  c(:, 6); mask = p<0.05;
    pairs =  num2cell(c(mask, 1:2),2);
    sigstar(pairs, p(mask)); hold off

    pvVal = orderpvalues(c); 
    pvVal(pvVal > 0.05) = NaN;
    plotheatmap(categories,pvVal)

% L/D Ratio
Mean = grpstats(DiodeLDRatio, CountLabels, 'mean');
StdDev = grpstats(DiodeLDRatio, CountLabels, 'std');
Err = StdDev./sqrt(LengthVector);

figure(Theme='light'); %set(gcf, 'Units', 'inches', 'Position', [1, 1, 10, 10]);
    b = bar(Mean); hold on; % excluding pulse train because the fit was so bad
    b.FaceColor = 'flat';         
    for i = 1:numel(Mean), b.CData(i,:) = Colors(i,:); end
    errorbar(1:numel(Mean), Mean, Err, 'k.', 'LineWidth', 1.5);
    title('Relative Diode Damage Assay'),   ylabel('\DeltaDead/Live Cells (%)');
    set(gca, 'XTick', 1:length(categories), 'XTickLabels', categories);
    set(gca,'Fontsize', 24)

    [~,~,statsVal] = anova1(DiodeLDRatio, CountLabels,'off'); 
    c = multcompare(statsVal,'CriticalValueType','bonferroni','Display','off');
    p =  c(:, 6); mask = p<0.05;
    pairs =  num2cell(c(mask, 1:2),2);
    sigstar(pairs, p(mask)); hold off

  
    PW(1).PIdF = (Exp(1).PostPIintensity - Exp(1).PrePIintensity)./...
                  Exp(1).PrePIintensity;
    PW(1).RE = Exp(1).RE;
    
    PW(2).PIdF = (Exp(2).PostPIintensity - Exp(2).PrePIintensity)./...
                  Exp(2).PrePIintensity;
    PW(2).RE = Exp(2).RE;
    
    PW(3).PIdF = (vertcat(Exp(3).PostPIintensity,Exp(5).PostPIintensity) - ...
            vertcat(Exp(3).PrePIintensity,Exp(5).PrePIintensity))./...
            vertcat(Exp(3).PrePIintensity,Exp(5).PrePIintensity);
    PW(3).RE = vertcat(Exp(3).RE,Exp(5).RE);
    
    PW(4).PIdF = (vertcat(Exp(4).PostPIintensity,Exp(6).PostPIintensity) - ...
            vertcat(Exp(4).PrePIintensity,Exp(6).PrePIintensity))./...
            vertcat(Exp(4).PrePIintensity,Exp(6).PrePIintensity);
    PW(4).RE = vertcat(Exp(4).RE,Exp(6).RE);
    
    Titles = {'Control','350 \mus','2 ms','8 ms'};
    
     for i = 1:length(PW)
         PW(i).Mask = PW(i).PIdF > max(PW(1).PIdF);
     end

%% Probit Plots
RErange = 0:.01:20;         ActivationProb = 0:0.01:.9;
figure(Theme='light');          set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
for i = 3:4
    binsize = iqr(PW(i).RE)*2/nthroot(length(PW(i).RE),3); % Freedman-Diaconis
    mask = PW(i).Mask; %sum(mask);

    % Probit Creation
    [PW(i).Probit, PW(i).bins, PW(i).cdf, PW(i).Probitgof]...
        = ProbitAnalysis(PW(i).RE,find(mask),binsize,0);
    PW(i).ProbitCurve = feval(PW(i).Probit,RErange);
    PW(i).ProbitCI = predint(PW(i).Probit, RErange, 0.95,'functional','off');
    Probit = [PW(i).ProbitCurve PW(i).ProbitCI];
    conf = confint(PW(i).Probit, 0.95); % 95% confidence intervals 
    SE = (conf(2, 1)-conf(1, 1))/2; % standard error of b
    SE = [0 ,SE, -SE];
    H50(i,:) = PW(i).Probit.a+SE;

    mpeak(i,:)  = max(diff(Probit,1,1),[],1)/0.01;
    
    % Group Plotting
    plot(RErange,PW(i).ProbitCurve,LineWidth=3,Color=Colors(i,:)); hold on 
    scatter(PW(i).bins,PW(i).cdf,50,'filled','MarkerEdgeColor','k','MarkerFaceColor',Colors(i,:),'HandleVisibility','off')

end

xlabel('Radiant Exposure (J/cm^2)'), ylabel('Damage Probability')
legend(Titles(3:4),"Location","northwest","Color",'none','EdgeColor','none')
ylim([-0.05 1.1]),  
ax = gca; ax.FontSize = 24;
labels = repmat(1:2, [3 1])';

% H50
figure(Theme='light');
    set(gcf, 'Units', 'inches', 'Position', [1, 1, 7,7]);
    b = bar(H50(3:4,1)); hold on;
    b.FaceColor = 'flat';
    for i = 3:4, b.CData(i-2,:) = Colors(i,:); end    
    errorbar(1:2, H50(3:4,1), H50(3:4,2)-H50(3:4,1),  H50(3:4,3)-H50(3:4,1), 'k.', 'LineWidth', 1.5);
    set(gca, 'XTick', 1:length(Titles(3:4)), 'XTickLabels', Titles(3:4));
    title("50% Damage Probability"), ylabel('H_5_0 (J/cm^2)'); ylim([0 15])
    ax = gca; ax.FontSize = 24; ax.YGrid = "on";

    [~,p] = ttest2(H50(3,:),H50(4,:));
    pairs =  num2cell(1:2,2);
    sigstar(pairs, p); hold off

% Probit Slopes
figure(Theme='light');
    set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
    b = bar(mpeak(3:4,1)); hold on;
    b.FaceColor = 'flat';
    for i = 3:4, b.CData(i-2,:) = Colors(i,:); end    
    errorbar(1:2, mpeak(3:4,1), mpeak(3:4,3)-mpeak(3:4,1),  mpeak(3:4,2)-mpeak(3:4,1), 'k.', 'LineWidth', 1.5);
    title("Slope of Probit Curves"), ylabel('m_p_e_a_k (J/cm^2)^-^1');
    set(gca, 'XTick', 1:length(Titles(3:4)), 'XTickLabels', Titles(3:4));
    ax = gca; ax.FontSize = 24;

    [~,p] = ttest2(mpeak(3,:),mpeak(4,:),"Tail","both");
    pairs =  num2cell(1:2,2);
    sigstar(pairs, p); hold off

%% TLA vs Theoretical spot comparisons
meanH50 = reshape(H50comp(:,1),[2,2])';
H50err = reshape(H50comp(:,2),[2,2])'-meanH50;
figure(Theme='light')
set(gcf, 'Units', 'inches', 'Position', [2, 2, 10, 7]);
    b = bar(meanH50); colororder(Colors(3:4,:)), hold on
    errorbar([b.XEndPoints], meanH50(:), H50err(:), 'k.', 'LineWidth', 1.5);
    set(gca, 'XTickLabels', {'Theoretical','TLA'});
    title("Threshold Comparison"), ylabel('H_5_0 (J/cm^2)'); ylim([0 15])
    ax = gca; ax.FontSize = 24;
    groupLabels = repmat(sort([b.XEndPoints]), [numel(H50comp(1,:)) 1])';
    
    % Intragroup comps
    for i = 1:2
        k = i+(i-1);
        mask = [k k+1];
        [~,p] = ttest2(H50comp(mask(1),:),H50comp(mask(2),:));
        pairs =  num2cell(groupLabels(mask(1):mask(2)),2);
        sigstar(pairs, p); 
    end
    
    % Intragroup comps
    for i = 1:2
        mask = [i i+2];
        [~,p] = ttest(H50comp(mask(1),:),H50comp(mask(2),:));
        pairs =  num2cell(groupLabels([mask(1),mask(2)]),2);
        sigstar(pairs, p); 
    end
    legend(Titles(3:4),"Location","northeast","Color",'none','EdgeColor','none')
%%
function pmat = orderpvalues(compoutput)

    pmat = ones(max(compoutput(:,2)),max(compoutput(:,2))) * 2;
    
    for i = 1:size(compoutput,1)
        
        pmat(compoutput(i,1),compoutput(i,2)) = compoutput(i,6);
        pmat(compoutput(i,2),compoutput(i,1)) = compoutput(i,6);
        
    end
    pmat(pmat>1)=NaN;

end