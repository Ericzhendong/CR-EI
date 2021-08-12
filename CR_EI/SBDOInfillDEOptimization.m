function [poptm foptm] = SBDOInfillDEOptimization(infillDEObj, objOPTS, ...
    normSpace, popsize, NbVariables)

% function [poptm foptm] = srgtsSBDOInfillDEOptimization(infillDEObj, objOPTS, ...
%     normSpace, popsize, NbVariables)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run differential evolution (optimizing the infill criterion)

% itermax = 5; % maximum of function evaluations = itermax*popsize
% foptm = Inf;
% for c1 = 1 : 1
itermax = 200; % maximum of function evaluations = itermax*popsize
foptm = Inf;



% for c1 = 1 : 4
%     [ptemp ftemp] = OPTMDE(infillDEObj, -Inf, NbVariables, ...
%         normSpace(1,:), normSpace(2,:), objOPTS, popsize, itermax, ...
%         0.8, 0.8, 7, 0); % DE parameters
%     if ftemp < foptm
%         poptm = ptemp;
%     end
% end

% [ptemp ftemp] = OPTMDE(infillDEObj, -Inf, NbVariables, ...
%     normSpace(1,:), normSpace(2,:), objOPTS, popsize, itermax, ...
%     0.8, 0.8, 7, 0); % DE parameters

% function [bestmem, bestval, nfeval] = DE(fname, D, XVmin, XVmax, NP, itermax)
[poptm foptm] = DE(infillDEObj,NbVariables,normSpace(1,:),normSpace(2,:),popsize,itermax);


% [poptm ftemp] = OPTMDE(infillDEObj, -Inf, NbVariables, ...
%         normSpace(1,:), normSpace(2,:), objOPTS, popsize, itermax, ...
%         0.8, 0.8, 7, 0); % DE parameters


return
