function [Area, radius] = SpotSizeCalc(NA,n1,n2,FiberDiameter,FiberDistance,CoverSlipThickness)
% all distances in um, n1 should be 1 for air and 1.33 for water, n2 should
% be 1.51 for glass
    FiberDiameter = FiberDiameter*(1-1/exp(2));
    FiberDivergenceAngle = asin(NA/n1); % divergence angle in air (n = 1)
    GlassDivergenceAngle = asin(sin(FiberDivergenceAngle)/n2); % snells law for continuing divergence in glass
    halfwidth = FiberDiameter/2 + ...
                FiberDistance*tan(FiberDivergenceAngle) + ...
                CoverSlipThickness*tan(GlassDivergenceAngle); % approximate beam half width in um
end