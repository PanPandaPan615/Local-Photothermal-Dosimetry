function  Data = AppendPWData(Dish)
    idx = (1:length(Dish(1).MeanF)-1);

    Data.PulseEnergy =    [];    Data.RE =     [];    Data.Energy = [];
    Data.Area =   [];

    Data.MeanF =   [];    Data.Std =     [];    Data.MaxF =    [];
    Data.Entropy = [];    Data.Skw =     [];    Data.Uni =     [];

    Data.Centroid = [];    Data.Rho =       [];    Data.Theta =     [];

    for i = 1:length(Dish)
        Data.PulseEnergy = [Data.PulseEnergy Dish(i).PulseEnergy];
        Data.RE =          [Data.RE Dish(i).RE];
        Data.Energy =      [Data.Energy Dish(i).Energy];
        Data.Area =        [Data.Area Dish(i).Area];
    
        Data.MeanF =   [Data.MeanF Dish(i).MeanF(idx,:)];
        Data.Std =     [Data.Std Dish(i).Std(idx,:)];
        Data.MaxF =    [Data.MaxF Dish(i).MaxF(idx,:)];
        Data.Entropy = [Data.Entropy Dish(i).Entropy(idx,:)];
        Data.Skw =     [Data.Skw Dish(i).Skw(idx,:)];
        Data.Uni =     [Data.Uni Dish(i).Uni(idx,:)];
    
        Data.Centroid = [Data.Centroid Dish(i).Centroid];
        Data.Rho =       [Data.Rho   Dish(i).Rho];
        Data.Theta =     [Data.Theta Dish(i).Theta];

    end

end