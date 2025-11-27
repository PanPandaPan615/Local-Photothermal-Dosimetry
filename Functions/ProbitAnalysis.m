function [Probit, bins, cdf, gof] = ProbitAnalysis(RadiantEnergy,Activated,binsize,plotit)
    % Radiant Energy
    bins  = 0:binsize:(round(max(RadiantEnergy),1)+binsize);
    binnedRE = discretize(RadiantEnergy,bins);
    binnedPass = discretize(RadiantEnergy(Activated),bins);
    bins = bins(1:end-1);
    
    for i = 1:length(bins)
        trials(i) = sum(binnedRE == i);
        pass(i) = sum(binnedPass == i);
    end
    
    cdf = pass./trials;
    nanidx = find(isnan(cdf));
    cdf(nanidx) = [];    bins(nanidx) = [];
    %ErrorFnt = fittype('d + (a - d) ./ (1 + exp(-b * (x - c)))');

    ErrorFnt = '0.5*(1+erf((x-a)/(b*sqrt(2))))';
    [Probit, gof] = fit(bins',cdf',ErrorFnt, 'StartPoint',[mean(bins), 1]);
 
end