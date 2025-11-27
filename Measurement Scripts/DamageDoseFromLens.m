clc, clear, close all
addpath(genpath(cd))
selpath = SelectFolders;

tic
for ExpFolders = 1:length(selpath)

    cd(selpath{ExpFolders});

    T = readtable("Experiment Settings.xlsx");
    ExpParameters = table2array(T);
    
    DishFolders = dir('Dish*');
    DishFolders = DishFolders([DishFolders.isdir]);
    
    res = 0.288798; % Optical Resolution um/px
    hwx = 241;%measured halfwidth of spotsize in um
    %hwx = SpotSizeCalc(0.22,1,1.51,400,200,100); % theoretical halfwidth    
    hwy = hwx;
    PW = [0, 0.35, 0.35, 2, 2, 2, 8, 8, 8];
    
    [X,Y] = meshgrid((0:2303)*res); % in um
    data = zeros(size(X,1),size(Y,2),2);
    data(:,:,1) = X; data(:,:,2) = Y;
    
    for D = 1:length(DishFolders)
        % Get IR fiducial Information
            cd("Fiducials\")
                Images = dir([pwd '/*.tif']); % load names of all .tif files in directory in struct
                FitInfo = ThermalLensFit(Images(D).name,res);
            cd ..
        
        % Create Maps
            PulseE = ExpParameters(:,end);
            nanidx = find(isnan(PulseE));
            PulseE(nanidx) = [];
            ux = FitInfo.meanx;  uy = FitInfo.meany;
            %hwx = FitInfo.hwx;  hwy = FitInfo.hwy; % optional halfwidth from lens
            H0 =  2*PulseE./(pi*hwx*hwy);
            E0 = H0*res^2; % mJ/um^2*um^2 = mJ
            DC = 0;
            
            EnergyModel = zeros(size(X,1),size(Y,2),length(PulseE));
            
                for i = 1:length(PulseE)
                    x0 = [E0(i),ux,hwx,uy,hwy, DC]; % gauss parameters [Amplitude, centerx, 1/e^2 x, centery, 1/e^2 y]
                    EnergyModel(:,:,i) = Gauss2D(x0,data); % Fitted Gauss
                end
       
        % Measure Energy Dose at cell soma        
            cd(DishFolders(D).name);
            cd("MergedBinary\");
                Binaries = dir('Binaries*');%   Load binary masks
                [EnergyTemp, ROI] = MeasureExposure(Binaries, EnergyModel, FitInfo, res);
            cd ..
        
        % Organize Data
            CountData = readtable("MergedCount.csv");
            IntensityData = readtable("MergedIntensity.csv");
            DishData = OrganizeData(EnergyModel, EnergyTemp, ROI, PulseE,...
                                    PW, CountData, IntensityData, FitInfo);
            cd ..
        
        % Save the data
        filename = DishFolders(D).name;
        filename(strfind(DishFolders(D).name,' ')) = [];
        filename = strcat(filename,'.mat');
        
        save(filename,"DishData");
        
    end
end
toc

%%
function paths = SelectFolders
    i = 1; folders = 1;
    while folders ~=0
        folders = uigetdir;
        if folders ~=0
            paths{i} = folders;
        end
        i = i+1;
    end

end

function FitInfo = ThermalLensFit(Image,res)

        I = double(imread(Image)); % I is image
        
        rowprof = mean(I,1)'; colprof = mean(I,2); % rough profile of gaussian
        
        fr = fit(res*(1:length(I))', rowprof,'gauss1'); % mean row profile for guess
        fc = fit(res*(1:length(I))', colprof,'gauss1'); % mean column profile for guess
        
        x0 = [(fr.a1+fc.a1),fr.b1,fr.c1/2,fc.b1,fc.c1/2, 0]; % initial guesses [Amplitude, centerx, 1/e^2 x, centery, 1/e^2 x, DC offset]
        [X,Y] = meshgrid((1:size(I))*res);
        data = zeros(size(X,1),size(Y,2),2);
        data(:,:,1) = X; data(:,:,2) = Y;
        [x,~,~,~] = lsqcurvefit(@Gauss2D,x0,data,I); % Least Squares Fit [x,~,residual,exitflag]
        %xerr = abs(x-x0)./x*100; % Error in initial guesses
        
        %store parameters in data structure
        FitInfo.name = Image; % 
        FitInfo.amplitude = x(1); % somewhat irrelavant thing to save
        FitInfo.meanx = x(2); % location of center in x
        FitInfo.hwx = x(3); %  Halfwidth in x
        FitInfo.meany = x(4); % location  of center in y
        FitInfo.hwy = x(5);% Halfwidth in y
        FitInfo.DC = x(6);% Halfwidth in y
        
        F = Gauss2D(x,data); % Fitted Gauss
    
        %imwrite(uint16(F),strcat('Fit_',Images(n).name),'tiff') % write fit to .tif
        
        FitInfo.Rsqr = diag(corrcoef(I,F),-1)^2; % calculate R^2 of fit 
 
    %writetable(struct2table(FitInfo),'FitInformation.csv'); % store fit information in table
end

function F = Gauss2D(x,data)
 F = x(1)*exp( -( 2*(data(:,:,1)-x(2)).^2/(x(3)^2) + 2*(data(:,:,2)-x(4)).^2/(x(5)^2) ) + x(end));
end

function halfwidth = SScalc(NA,n1,n2,FiberDiameter,FiberDistance,CoverSlipThickness)
% all distances in um, n1 should be 1 for air and 1.33 for water, n2 should
% be 1.51 for glass

    FiberDivergenceAngle = asin(NA/n1); % divergence angle in air (n = 1)
    GlassDivergenceAngle = asin(sin(FiberDivergenceAngle)/n2); % snells law for continuing divergence in glass
    halfwidth = FiberDiameter/2 + ...
                FiberDistance*tan(FiberDivergenceAngle) + ...
                CoverSlipThickness*tan(GlassDivergenceAngle); % approximate beam half width in um
end

function [EnergyMeasurment, ROI] = MeasureExposure(Binaries, EnergyMaps,FitInfo,res)
ROI.Binary = zeros(size(EnergyMaps));

    for n = 1:length(Binaries)    
        
        BinTemp = imread(Binaries(n).name);
        B = imbinarize(BinTemp);

        ROI.Binary(:,:,n) = B;
        I = EnergyMaps(:,:,n);
        ROIs = bwlabel(B',8)';% Label ROIs

        props = regionprops(ROIs, I, 'WeightedCentroid','Area','PixelValues');
        
        Area = cat(1, props.Area)*(res^2); % Area is returned at number of pixels * res^2 um^2/px
        Energy = zeros(1,length(Area))';
        for i = 1:length(Area)
            Energy(i) = sum(props(i).PixelValues); % mJ
        end
        
        EnergyMeasurment(n).RadiantEnergy = Energy./Area*1e5;
        EnergyMeasurment(n).Area = Area;
        EnergyMeasurment(n).Energy = Energy;

        ROIdetails(n).centroids = cat(1, props.WeightedCentroid);
        muIR = [FitInfo.meanx/res FitInfo.meany/res]; % pixel location
        if isempty(ROIdetails(n).centroids)
            adjusted_centroid = muIR;
        else
            adjusted_centroid = muIR - ROIdetails(n).centroids;
        end
        ROIdetails(n).rho = sqrt(adjusted_centroid(:,1).^2 + adjusted_centroid(:,2).^2);
        ROIdetails(n).theta = atan2(adjusted_centroid(:,2),adjusted_centroid(:,1));

    end

    ROI.CentroidData = ROIdetails;
end

function [DishData] = OrganizeData(EnergyModel, EnergyTemp, ROI, PulseEnergy,...
                                    PW, CountData, IntensityData, FitInfo)

for n = 1:max(CountData{:,"MultiPointIndex"})

    MPidx = IntensityData{:,"MultiPointIndex"} == n;
    RawData = IntensityData(MPidx,:);

    Measurements(n).PrePIintensity  = RawData{:,"PreExposure_PI_Intensity"};
    Measurements(n).PostPIintensity = RawData{:,"PostExposure_PI_Intensity2"};
    Measurements(n).PrePICount      = CountData{n,"PreExposure_PI"};
    Measurements(n).CalcienCount    = CountData{n,"Calcien"};
    Measurements(n).PostPICount     = CountData{n,"PostExposure_PI"};

    Measurements(n).RadiantEnergy = EnergyTemp(n).RadiantEnergy;
    Measurements(n).Area = EnergyTemp(n).Area;
    Measurements(n).Energy = EnergyTemp(n).Energy;
    Measurements(n).PulseEnergy = PulseEnergy(n);
    Measurements(n).PW = PW(n);

end

    DishData.Measurements = Measurements;

    DishData.ROIdata = ROI;
    DishData.EnergyModel = EnergyModel;
    DishData.FitInfo = FitInfo;
    
end