function [r,f] = regain(f,E,toi,blperiod)    
    r         = corr(f(toi)',E(1,toi)').^2;    
    f         = utils.baseline(f,blperiod,1);
end