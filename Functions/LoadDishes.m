function [DishID Dish GoF] = LoadDishes()
    DishID = dir('Dish*.mat');
    for i = 1:length(DishID)
       temp(i)= load(DishID(i).name);
       Dish(i) = temp(i).DishData;
    
       GoF(i,1) = Dish(i).FitInfo.Rsqr;
       GoF(i,2) = Dish(i).FitInfo.hwx;
       GoF(i,3) = Dish(i).FitInfo.hwy;
    end
end