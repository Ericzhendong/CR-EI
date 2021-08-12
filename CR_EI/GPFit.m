function [FIT] = GPFit(OPT)

%% Initial Settings

X = OPT.X;
Y = OPT.Y;
[NbPoints NbVariables] = size(X);
likfunc = @likGauss;
covfunc = OPT.covfunc;

%% Normalization

mx = mean(X);
sx = std(X);
my = mean(Y);
sy = std(Y);
nparx  = [mx;sx];
npary  = [my;sy];

X_Norm = (X-repmat(mx,size(X,1),1))./repmat(sx,size(X,1),1); 
Y_Norm = (Y-my)/sy;   % Fitting purpose

%% One-short hyper-parameter tuning
hyp0 = OPT.hyp;
hyp_Pre  = minimize(hyp0,@gp,-200,@infGaussLik,[],covfunc,likfunc,X_Norm,Y_Norm);
hyp = hyp_Pre;

% XVmin = [0.5*log(0.1)*ones(NbVariables,1);log(0.5);log(sqrt(1.0e-9))];
% XVmax = [0.5*log(10)*ones(NbVariables,1);log(1.5);log(sqrt(1.0e-6))];
% 
% % XVmin = [log(sqrt(1.0e-9));0.5*log(0.1)*ones(NbVariables,1);log(0.5)];
% % XVmax = [log(sqrt(1.0e-1));0.5*log(10)*ones(NbVariables,1);log(1.5)];
% 
% hyp = boxmin(hyp_Pre, @gp, XVmin, XVmax, @infGaussLik,[],covfunc,likfunc,X_Norm,Y_Norm);


% hyp.cov = [-0.9987 -0.0294]';
% hyp.lik = -13.7861;

%%
FIT.X_Norm = X_Norm;
FIT.Y_Norm = Y_Norm; 
FIT.nparx  = nparx;
FIT.npary  = npary;

FIT.NbPoints    = NbPoints;
FIT.NbVariables = NbVariables;

FIT.covfunc = covfunc;
FIT.hyp     = hyp;
FIT.likfunc = likfunc;

