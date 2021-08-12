% This document includes two parts:
% 1. Algorithm Testing with Different Compared Algorithms such as
% CR-EI,EGO,PEI,GEI and E3I
% 2. Results Analysis such as convergence history and boxplot
% 3. Note that testing of the baseline E3I can be very time-consuming, as it involes the
% process of Thompson sampling

%% Hartman4 Function

myFN    = @Hartman4;

designspace = [0 0 0 0;  % lower bound
               1 1 1 1]; % upper bound
          
ndv = size(designspace, 2);

maxNbCycles  = 40*size(designspace,2);
npoints = 6*size(designspace,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Testing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing 

for i = 1:20
    
    TestNumber = i;
    
    FuncName = 'Hartman4\';
    FilePath_mkdir = ['TestingResults\',FuncName];
    
    mkdir(FilePath_mkdir,num2str(TestNumber));
    
    
    %% Design of Experiments
    dLHS = lhsdesign(npoints,ndv,'criterion','maximin','iteration',1000);
    X = srgtsScaleVariable(dLHS, [zeros(1, ndv); ones(1, ndv)], ...
        designspace);
    Y = feval(myFN, X);

    DataSamples = [X Y];
    
    FilePath = ['TestingResults\',FuncName, num2str(TestNumber),'\DataSamples.txt'];
    save(FilePath,'DataSamples','-ascii');

    
    %% Fitting of GP Surrogate
    SOGP_OPT = SetOptions('GP',X, Y); 
    SOGP_FIT = GPFit(SOGP_OPT);

    %% Compared Algorithms
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CR-EI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x, y, stateCR_EI] = SBDODriver_CR_EI(myFN,designspace, designspace,... 
        SOGP_OPT, SOGP_FIT, [], @ExpectedImprovement, maxNbCycles, 'ExpectedImprovement',1);

    format long;
    Results_CR_EI = [stateCR_EI.Ypbs',stateCR_EI.neval];

    FilePath = ['TestingResults\',FuncName, num2str(TestNumber),'\Results_CR_EI.txt'];
    save(FilePath,'Results_CR_EI','-ascii');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EGO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [x, y, stateEGO] = SBDODriver_EGO(myFN,designspace, designspace,... 
        SOGP_OPT, SOGP_FIT, [], @ExpectedImprovement, maxNbCycles, 'ExpectedImprovement',1);

    format long;
    ResultsEGO = [stateEGO.Ypbs',stateEGO.neval];

    FilePath = ['TestingResults\',FuncName, num2str(TestNumber),'\ResultsEGO.txt'];
    save(FilePath,'ResultsEGO','-ascii');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PEI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x, y, statePEI] = SBDODriver_PEI(myFN,designspace, designspace,... 
        SOGP_OPT, SOGP_FIT, [], @ExpectedImprovement, maxNbCycles/2, 'ExpectedImprovement',2);

    format long;
    ResultsPEI = [statePEI.Ypbs',statePEI.neval];

    FilePath = ['TestingResults\',FuncName, num2str(TestNumber),'\ResultsPEI.txt'];
    save(FilePath,'ResultsPEI','-ascii');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GEI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x, y, stateGEI] = SBDODriver_GEI(myFN,designspace, designspace,... 
        SOGP_OPT, SOGP_FIT, [], @Generalized_ExpectedImprovement, maxNbCycles, 'Generalized_ExpectedImprovement',2);

    format long;
    ResultsGEI = [stateGEI.Ypbs',stateGEI.neval];

    FilePath = ['TestingResults\',FuncName, num2str(TestNumber),'\ResultsGEI.txt'];
    save(FilePath,'ResultsGEI','-ascii');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E3I%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    npointspercycle = 1;
    [x, y, stateE3I] = SBDODriver_E3I(myFN, designspace, designspace,... 
        SOGP_OPT, SOGP_FIT, [], @E3I, maxNbCycles, 'E3I',npointspercycle);

    format long;
    ResultsE3I = [stateE3I.Ypbs',stateE3I.neval];

    FilePath = ['TestingResults\',FuncName, num2str(TestNumber),'\ResultsE3I.txt'];
    save(FilePath,'ResultsE3I','-ascii');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Results Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results Analysis

Func_name = 'TestingResults\Hartman4\'

n = 20;
Ybest = -3.134500;

OptCR_EI = zeros(n,1);
OptEGO = zeros(n,1);
OptPEI = zeros(n,1);
OptGEI = zeros(n,1);
OptE3I = zeros(n,1);

ResultsCR_EI = zeros(maxNbCycles,n);
ResultsEGO = zeros(maxNbCycles,n);
ResultsPEI = zeros(maxNbCycles/2,n);
ResultsGEI = zeros(maxNbCycles,n);
ResultsE3I = zeros(maxNbCycles,n);

for i = 1:n
    
    TestNumber = i;
    
    % Load data
    
    FilePathCR_EI = [Func_name,num2str(TestNumber),'\Results_CR_EI.txt'];
    ResultsCR_EI_tmp  = load(FilePathCR_EI);
    
    FilePathEGO = [Func_name,num2str(TestNumber),'\ResultsEGO.txt'];
    ResultsEGO_tmp  = load(FilePathEGO);
    
    FilePathPEI = [Func_name,num2str(TestNumber),'\ResultsPEI.txt'];
    ResultsPEI_tmp  = load(FilePathPEI);
    
    FilePathGEI = [Func_name,num2str(TestNumber),'\ResultsGEI.txt'];
    ResultsGEI_tmp  = load(FilePathGEI);
    
    FilePathE3I = [Func_name,num2str(TestNumber),'\ResultsE3I.txt'];
    ResultsE3I_tmp  = load(FilePathE3I);
    
    % Best solution using boxplot
    OptCR_EI(i) = min(ResultsCR_EI_tmp(:,1)) - Ybest;
    OptEGO(i) = min(ResultsEGO_tmp(:,1)) - Ybest;
    OptPEI(i) = min(ResultsPEI_tmp(:,1)) - Ybest;
    OptGEI(i) = min(ResultsGEI_tmp(:,1)) - Ybest;
    OptE3I(i) = min(ResultsE3I_tmp(:,1)) - Ybest;

    % Convergence history
    ResultsCR_EI(:,i) = ResultsCR_EI_tmp(:,1)- Ybest;
    ResultsEGO(:,i) = ResultsEGO_tmp(:,1)- Ybest;   
    ResultsPEI(:,i) = ResultsPEI_tmp(:,1)- Ybest; 
    ResultsGEI(:,i) = ResultsGEI_tmp(:,1)- Ybest; 
    ResultsE3I(:,i) = ResultsE3I_tmp(:,1)- Ybest;   

    IterCR_EI1(:,i) = ResultsCR_EI_tmp(:,2);
    IterEGO1(:,i) = ResultsEGO_tmp(:,2);
    IterPEI1(:,i) = ResultsPEI_tmp(:,2);
    IterGEI1(:,i) = ResultsGEI_tmp(:,2);
    IterE3I1(:,i) = ResultsE3I_tmp(:,2);
     
end


%% Convergence history

ResultsCR_EI = log10(ResultsCR_EI);
ResultsEGO = log10(ResultsEGO);
ResultsPEI = log10(ResultsPEI);
ResultsGEI = log10(ResultsGEI);
ResultsE3I = log10(ResultsE3I);

MeanCR_EI = mean(ResultsCR_EI,2);
MeanEGO = mean(ResultsEGO,2);
MeanPEI = mean(ResultsPEI,2);
MeanGEI = mean(ResultsGEI,2);
MeanE3I = mean(ResultsE3I,2);

StdCR_EI = (std(ResultsCR_EI,0,2));
StdEGO = (std(ResultsEGO,0,2));
StdPEI = (std(ResultsPEI,0,2));
StdGEI = (std(ResultsGEI,0,2));
StdE3I = (std(ResultsE3I,0,2));

IterCR_EI = ceil(mean(IterCR_EI1,2));
IterEGO = ceil(mean(IterEGO1,2));
IterPEI = ceil(mean(IterPEI1,2));
IterGEI = ceil(mean(IterGEI1,2));
IterE3I = ceil(mean(IterE3I1,2));

logyCR_EI = MeanCR_EI;
logyEGO = MeanEGO;
logyPEI = MeanPEI;
logyGEI = MeanGEI;
logyE3I = MeanE3I;

figure();
f1 = plot(IterEGO,logyEGO,'-b*');
hold on;
f2 = plot(IterPEI,logyPEI,'-rp');
hold on;
f3 = plot(IterGEI,logyGEI,'-gs');
hold on;
f4 = plot(IterE3I,logyE3I,'-c+');
hold on;
f5 = plot(IterCR_EI,logyCR_EI,'-ko','MarkerFaceColor','k');
hold on;



ciub_EGO = MeanEGO(1:end) + StdEGO(1:end); 
cilb_EGO = MeanEGO(1:end) - StdEGO(1:end); 
xShaded = [IterEGO(1:end);sort(IterEGO(1:end),'descend')]; 
yShaded = [cilb_EGO;ciub_EGO(end:-1:1)];
patch(xShaded,yShaded,[.24 .73 .89],'EdgeColor','none','FaceColor','blue','FaceAlpha',0.1)


ciub_PEI = MeanPEI(1:end) + StdPEI(1:end); 
cilb_PEI = MeanPEI(1:end) - StdPEI(1:end); 
xShaded = [IterPEI(1:end);sort(IterPEI(1:end),'descend')]; 
yShaded = [cilb_PEI;ciub_PEI(end:-1:1)];
patch(xShaded,yShaded,[.24 .73 .89],'EdgeColor','none','FaceColor','red','FaceAlpha',.1);


ciub_GEI = MeanGEI(1:end) + StdGEI(1:end); 
cilb_GEI = MeanGEI(1:end) - StdGEI(1:end); 
xShaded = [IterGEI(1:end);sort(IterGEI(1:end),'descend')]; 
yShaded = [cilb_GEI;ciub_GEI(end:-1:1)];
patch(xShaded,yShaded,[.24 .73 .89],'EdgeColor','none','FaceColor','green','FaceAlpha',.1);


ciub_E3I = MeanE3I(1:end) + StdE3I(1:end); 
cilb_E3I = MeanE3I(1:end) - StdE3I(1:end); 
xShaded = [IterE3I(1:end);sort(IterE3I(1:end),'descend')]; 
yShaded = [cilb_E3I;ciub_E3I(end:-1:1)];
patch(xShaded,yShaded,[.24 .73 .89],'EdgeColor','none','FaceColor','cyan','FaceAlpha',.1);


ciub_CR_EI = MeanCR_EI(1:end) + StdCR_EI(1:end); 
cilb_CR_EI = MeanCR_EI(1:end) - StdCR_EI(1:end); 
xShaded = [IterCR_EI(1:end);sort(IterCR_EI(1:end),'descend')]; 
yShaded = [cilb_CR_EI;ciub_CR_EI(end:-1:1)];
patch(xShaded,yShaded,[.24 .73 .89],'EdgeColor','none','FaceColor','black','FaceAlpha',.2);


legend([f1,f2,f3,f4, f5],'EGO','PEI','GEI','E3I','CR-EI');
legend('boxoff');

xlabel('Number of function evaluations');
ylabel('Log10 R_t');
set(gca,'fontsize',10,'fontweight','bold');



%% box plot
figure();
C1 = zeros(size(OptEGO));
C2 = ones(size(OptPEI));
C3 = 2*ones(size(OptGEI));
C4 = 3*ones(size(OptE3I));
C5 = 4*ones(size(OptCR_EI));

G = [C1' C2' C3' C4' C5'];
X = log10([OptEGO' OptPEI' OptGEI' OptE3I', OptCR_EI']);

boxplot(X,G,'Labels',{'EGO','PEI','GEI','E3I','CR-EI'});
set(gca,'fontsize',10,'fontweight','bold');

