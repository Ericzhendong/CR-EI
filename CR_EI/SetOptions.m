function OPT = SetOptions(Flag,X, Y, nTasks,nIndex)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options

switch Flag
    
    case 'GP'
        
        OPT.SRGT = Flag;
        [npoints nvariables] = size(X);  
        
        OPT.X = X;
        OPT.Y = Y;        
        OPT.covfunc = 'covSEard';
%         OPT.hyp.cov = [log(npoints^(-1/nvariables))*ones(nvariables, 1)/10; 0];
        OPT.hyp.cov = [0.5*log(0.5*ones(nvariables, 1));0.1]; 
%         OPT.hyp.cov = [0.5*ones(nvariables, 1);0.1]; % 20210105
        OPT.hyp.lik = log(sqrt(1.0e-3));

   
    case 'MTGP'
        
        OPT.SRGT = Flag;
        [npoints nvariables] = size(X); 

        OPT.nTasks  = nTasks;
        OPT.iTask   = 0;
        OPT.X       = X;
        OPT.Y       = Y; 
        
        OPT.nIndex    = nIndex;  
        
end

