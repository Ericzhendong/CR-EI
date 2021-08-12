function [stateEGO npoints ndv popsize] = SBDOMakeState(OPT)
% function [stateEGO npoints ndv popsize] = srgtsSBDOMakeState(srgtOPT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

stateEGO.Iteration    = 0;
stateEGO.ReasonToStop = [];


switch OPT.SRGT
    case 'GP'
        stateEGO.Ypbs = min(OPT.Y);
        [npoints ndv] = size(OPT.X);
        
    case 'MTGP'
        
        nTasks     = OPT.nTasks;
        nIndex     = OPT.nIndex;
        TransCodes = OPT.TransCodes;

        disp(['Iteration 0:']);
        
        for i = 1:nTasks
            
            Index_tmp = sum(nIndex(1:i-1));  
%             stateEGO.Ypbs(1,i) = min(OPT.Y(Index_tmp+1:Index_tmp+nIndex(i)));
            Ymin_tmp  = min(OPT.Y(Index_tmp+1:Index_tmp+nIndex(i)));
            %function y = Transformation(nstat,y,TransCode)
            stateEGO.Ypbs(1,i) = min(OPT.Y(Index_tmp+1:Index_tmp+nIndex(i)));
            disp(['Best function value for Task ',num2str(i),': ',num2str(stateEGO.Ypbs(1,i))]);

        end
               
        [npoints ndv] = size(OPT.X);  
        ndv = ndv - 1;  %The last term is index
    
    
    
    case 'Swersky'
        
        nTasks     = OPT.nTasks;
        nIndex     = OPT.nIndex;
        TransCodes = OPT.TransCodes;

        disp(['Iteration 0:']);
        
        for i = 1:nTasks
            
            Index_tmp = sum(nIndex(1:i-1));  
            Ymin_tmp  = min(OPT.Y(Index_tmp+1:Index_tmp+nIndex(i)));
            stateEGO.Ypbs(1,i) = Transformation(1,Ymin_tmp,TransCodes(i));
            disp(['Best function value for Task ',num2str(i),': ',num2str(stateEGO.Ypbs(1,i))]);

        end
               
        [npoints ndv] = size(OPT.X);  
        ndv = ndv - 1;  %The last term is index 
        
end


stateEGO.X = [];
stateEGO.Y = [];



if ndv < 10
    popsize = max(20*ndv, 50);
else
    popsize = 100;
end

return
