% This document intoduces how to run the CR-EI algorithm, which consists of
% the following four steps:
% 1. Initial Settings;
% 2. Design of Experiments
% 3. Fitting of GP surrogate
% 4. Algorithm Run 

myFN    = @Hartman3;

%% Initial Settings
FuncName = 'Hartman3\';
FilePath_mkdir = ['TestingResults\',FuncName];

mkdir(FilePath_mkdir);


designspace = [0 0 0;  % lower bound
               1 1 1]; % upper bound
          
ndv = size(designspace, 2);

maxNbCycles  = 40*size(designspace,2);
npoints = 6*size(designspace,2);


%% Design of Experiments
dLHS = lhsdesign(npoints,ndv,'criterion','maximin','iteration',1000);
X = srgtsScaleVariable(dLHS, [zeros(1, ndv); ones(1, ndv)], ...
    designspace);
Y = feval(myFN, X);

DataSamples = [X Y];


FilePath = ['TestingResults\',FuncName, '\DataSamples.txt'];
save(FilePath,'DataSamples','-ascii');


%% Fitting of GP Surrogate
SOGP_OPT = SetOptions('GP',X, Y); 
SOGP_FIT = GPFit(SOGP_OPT); 


%% Algorithm Run 
[x, y, stateCR_EI] = SBDODriver_CR_EI(myFN,designspace, designspace,... 
    SOGP_OPT, SOGP_FIT, [], @ExpectedImprovement, maxNbCycles, 'ExpectedImprovement',1);

format long;
Results_CR_EI = [stateCR_EI.Ypbs',stateCR_EI.neval];

FilePath = ['TestingResults\',FuncName, '\Results_CR_EI.txt'];
save(FilePath,'Results_CR_EI','-ascii');