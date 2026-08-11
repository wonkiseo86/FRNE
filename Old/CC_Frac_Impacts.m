function [Target0_grid, ref_recover, Proj_raw1, Proj_raw2, Proj_raw3, Proj_raw4, Proj_raw5,Sum_stat] = CC_Frac_Impacts(Target0_grid,Target0,Tef0_gw1,Wght)
T0 = 34; Temp = 0;           
for i =1:T0
    Temp0 = Target0(i,:)'/mean(Target0(i,:));
    Temp = Temp + Temp0;
end
ref = (1/T0)*Temp;
refclr = log(ref/geomean(ref)); ref_recover = exp(refclr)/mean(exp(refclr));

Proj_raw_stack = zeros(length(ref_recover),30);
for i = 1:30
    Temp_Proj_clr = refclr + Wght(i)*Tef0_gw1(:,1); Temp_Proj_raw = exp(Temp_Proj_clr)/mean(exp(Temp_Proj_clr)); 
    Proj_raw_stack(:,i) = Temp_Proj_raw;
end
Sum_stat = KS_Desc_stat(Target0_grid,[ref_recover'; Proj_raw_stack']);
Proj_raw1 = Proj_raw_stack(:,7); 
Proj_raw2 = Proj_raw_stack(:,13); 
Proj_raw3 = Proj_raw_stack(:,19);
Proj_raw4 = Proj_raw_stack(:,24); 
Proj_raw5 = Proj_raw_stack(:,30);
