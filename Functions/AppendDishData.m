function  Appended = AppendDishData(Dish,indicies)
    
Dish.Appended.PulseEnergy = []; Dish.Appended.RE =     [];
Dish.Appended.Energy = []; Dish.Appended.Area =   [];

Dish.Appended.MeanF =   []; Dish.Appended.Std =     [];
Dish.Appended.MaxF =    []; Dish.Appended.Entropy = [];
Dish.Appended.Skw =     []; Dish.Appended.Uni =     [];

Dish.Appended.Centroid =    []; Dish.Appended.Rho =         [];
Dish.Appended.Theta =       [];

    for i = indicies        
       PulseE = (ones(length(Dish.EnergyMeasurement.Measurement(i).RadiantEnergy),1)*Dish.EnergyMeasurement.PulseEnergy(i))';
    
       Dish.Appended.PulseEnergy = [Dish.Appended.PulseEnergy PulseE];
       Dish.Appended.RE =     [Dish.Appended.RE Dish.EnergyMeasurement.Measurement(i).RadiantEnergy'];
       Dish.Appended.Energy = [Dish.Appended.Energy Dish.EnergyMeasurement.Measurement(i).Energy'];
       Dish.Appended.Area =   [Dish.Appended.Area Dish.EnergyMeasurement.Measurement(i).Area'];
    
       Dish.Appended.MeanF =   [Dish.Appended.MeanF Dish.TimeMeasurement(i).meanintensity];
       Dish.Appended.Std =     [Dish.Appended.Std Dish.TimeMeasurement(i).stdintensity];
       Dish.Appended.MaxF =    [Dish.Appended.MaxF Dish.TimeMeasurement(i).stdintensity];
       Dish.Appended.Entropy = [Dish.Appended.Entropy Dish.TimeMeasurement(i).stdintensity];
       Dish.Appended.Skw =     [Dish.Appended.Skw  Dish.TimeMeasurement(i).cellskewness];
       Dish.Appended.Uni =     [Dish.Appended.Uni Dish.TimeMeasurement(i).celluniformity];
    
       Dish.Appended.Centroid =    [Dish.Appended.Centroid Dish.ROIdata.CentroidData(i).centroids'];
       Dish.Appended.Rho =         [Dish.Appended.Rho Dish.ROIdata.CentroidData(i).rho'];
       Dish.Appended.Theta =       [Dish.Appended.Theta Dish.ROIdata.CentroidData(i).theta'];
    end
    Appended = Dish.Appended;
end