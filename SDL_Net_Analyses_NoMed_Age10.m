function SDL_Net_Analyses_NoMed_Age10(SDL)

% Graph theory analyses for correlation and regression matrix
% Input
% ------ SDL, all important parameters
% Output
% ------ SDL, with newly added SDL.Out, the structure containing all results
%% add path of BCT toolbox
npath = fullfile(SDL.BoxPath,'BCT','2017_01_15_BCT'); % path of the BCT toolbox
addpath(genpath(npath));

%% Load data
fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
fn = fullfile(fdir,['WC_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: wiring cost & thresholds in %s\n',fn);
Thr = [W.Thr_abs;W.Thr_pos;W.Thr_neg]; % Thresholds
WC  = [W.WC_abs; W.WC_pos; W.WC_neg];  % wiring costs

fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Matrix in %s\n',fn);
% MRHO = rand(148*5,148,1000); % to test the procedure

if strcmp(SDL.ana_type{1},'corr') ||strcmp(SDL.ana_type{1},'partialcorr') % correlation matrix
    MMc = MRHO;
elseif strcmp(SDL.ana_type{1},'regr') || strcmp(SDL.ana_type{1},'regrcov') % regression matrix
    MMc = MB;
else
end
atype = 'und'; % connection without direction

%% Nodal Measures
% (1) Nodal Degree Centrality
[dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@degrees_und,MMc,WC,Thr);
SDL.Node_Degree = struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);
% (2) Nodal Betweenness Centrality
[dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@betweenness_bin,MMc,WC,Thr);
SDL.Node_Betweenness = struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);
% (3) Nodal Closeness Centrality (i.e. the reciprocal of distance, i.e. close_g1vg2 = 1/distance1-1/distance2)
[dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@distance_bin,MMc,WC,Thr);
SDL.Node_Closeness = struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);
% (4) Eigenvector Centrality
[dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@eigenvector_centrality_und,MMc,WC,Thr);
SDL.Node_Eigenvector = struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);

% Mc_abs
SDL.val_Mc_abs = [SDL.Node_Degree.val_abs(:,1),SDL.Node_Betweenness.val_abs(:,1),SDL.Node_Closeness.val_abs(:,1),SDL.Node_Eigenvector.val_abs(:,1)]; 
SDL.dif_Mc_abs = [SDL.Node_Degree.dif_abs(:,1),SDL.Node_Betweenness.dif_abs(:,1),SDL.Node_Closeness.dif_abs(:,1),SDL.Node_Eigenvector.dif_abs(:,1)];
SDL.p_Mc_abs   = [SDL.Node_Degree.p_abs(:,1),  SDL.Node_Betweenness.p_abs(:,1),  SDL.Node_Closeness.p_abs(:,1),  SDL.Node_Eigenvector.p_abs(:,1)];

% Mc_pos
SDL.val_Mc_pos = [SDL.Node_Degree.val_pos(:,1),SDL.Node_Betweenness.val_pos(:,1),SDL.Node_Closeness.val_pos(:,1),SDL.Node_Eigenvector.val_pos(:,1)]; 
SDL.dif_Mc_pos = [SDL.Node_Degree.dif_pos(:,1),SDL.Node_Betweenness.dif_pos(:,1),SDL.Node_Closeness.dif_pos(:,1),SDL.Node_Eigenvector.dif_pos(:,1)];
SDL.p_Mc_pos   = [SDL.Node_Degree.p_pos(:,1),  SDL.Node_Betweenness.p_pos(:,1),  SDL.Node_Closeness.p_pos(:,1),  SDL.Node_Eigenvector.p_pos(:,1)];

% Mc_neg
SDL.val_Mc_neg = [SDL.Node _Degree.val_neg(:,1),SDL.Node_Betweenness.val_neg(:,1),SDL.Node_Closeness.val_neg(:,1),SDL.Node_Eigenvector.val_neg(:,1)]; 
SDL.dif_Mc_neg = [SDL.Node_Degree.dif_neg(:,1),SDL.Node_Betweenness.dif_neg(:,1),SDL.Node_Closeness.dif_neg(:,1),SDL.Node_Eigenvector.dif_neg(:,1)];
SDL.p_Mc_neg   = [SDL.Node_Degree.p_neg(:,1),  SDL.Node_Betweenness.p_neg(:,1),  SDL.Node_Closeness.p_neg(:,1),  SDL.Node_Eigenvector.p_neg(:,1)];

% %% Edge Measures
% % (1) Edge Betweenness Centrality
% [dif_abs,p_abs,dif_pos,p_pos,dif_neg,p_neg] = SDL_Graph_output(@edge_betweenness_bin,MMc,WC,Thr);
% SDL.Edge_Betweenness = struct('dif_abs',dif_abs,'p_abs',p_abs,'dif_pos',dif_pos,'p_pos',p_pos,'dif_neg',dif_neg,'p_neg',p_neg);

%% save results
fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
fn = fullfile(fdir,['Reports_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
save(fn,'SDL'); fprintf('Saved: Reports in %s\n',fn);

% end
end

% Graph theory comparisons between M1 and M2 for abs, pos & neg mediations
function [dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(myfun,Mci,WCC,Thrr)
% Input
% ------ myfun, Matlab handle function, e.g. @degrees_und
% ------ Mci,  correlation or regression coefficients
% ------ WC & Thr, wiring cost and thresholds matrix
%                  odd  rows -- group 1
%                  even rows -- group 2
%                  rows 1-2, abs; rows 3-4, pos; rows 5-6, neg
% Output
% ------ dif, p & val, between-group differences values, p values and raw values per group for abs, pos and neg

N = 5; % number of age-ranges
k = size(Mci,2); % the big matrix Mci can be seperated into a few small matrix (k x k)
g11 = [0*k+1:1*k]; g12 = [1*k+1:2*k]; g13 = [2*k+1:3*k];  g14 = [3*k+1:4*k]; g18 = [4*k+1:5*k]; % index of group 1
g21 = [5*k+1:6*k]; g22 = [6*k+1:7*k]; g23 = [7*k+1:8*k];  g24 = [8*k+1:9*k]; g28 = [9*k+1:10*k]; % index of group 2

% absolute values
Thr = Thrr(1:N*2,:); WC = WCC(1:N*2,:);
% ThC  = Thr(find(WC==max(WC(:,1)))); % the threshold corresponding to the maximum value of minimum wring cost
% Mc  = double(abs(Mci)>ThC(1));
Mc = SDL_threshold_proportional(abs(Mci),max(WC(:,1)));
Mc = double(Mc>0);
[dif_abs(:,:,1),p_abs(:,:,1),val_abs(:,:,1),val_abs(:,:,2)] = SDL_Graph_Analyses(myfun,Mc(g11,:,:), Mc(g21,:,:));
[dif_abs(:,:,2),p_abs(:,:,2),val_abs(:,:,3),val_abs(:,:,4)] = SDL_Graph_Analyses(myfun,Mc(g12,:,:), Mc(g22,:,:));
[dif_abs(:,:,3),p_abs(:,:,3),val_abs(:,:,5),val_abs(:,:,6)] = SDL_Graph_Analyses(myfun,Mc(g13,:,:), Mc(g23,:,:));
[dif_abs(:,:,4),p_abs(:,:,4),val_abs(:,:,7),val_abs(:,:,8)] = SDL_Graph_Analyses(myfun,Mc(g14,:,:), Mc(g24,:,:));
%[dif_abs(:,:,5),p_abs(:,:,5),val_abs(:,:,9),val_abs(:,:,10)] = SDL_Graph_Analyses(myfun,Mc(g15,:,:), Mc(g25,:,:));
%[dif_abs(:,:,6),p_abs(:,:,6),val_abs(:,:,11),val_abs(:,:,12)] = SDL_Graph_Analyses(myfun,Mc(g16,:,:), Mc(g26,:,:));
%[dif_abs(:,:,7),p_abs(:,:,7),val_abs(:,:,13),val_abs(:,:,14)] = SDL_Graph_Analyses(myfun,Mc(g17,:,:), Mc(g27,:,:));
[dif_abs(:,:,5),p_abs(:,:,5),val_abs(:,:,9),val_abs(:,:,10)] = SDL_Graph_Analyses(myfun,Mc(g18,:,:), Mc(g28,:,:));

% pos values
Thr = Thrr(2*N+1:4*N,:); WC = WCC(2*N+1:4*N,:);
% ThC  = Thr(find(WC==max(WC(:,1)))); % the threshold corresponding to the maximum value of minimum wring cost
% Mc  = double(Mci>ThC(1));
Mc = SDL_threshold_proportional(Mci,max(WC(:,1)));
Mc = double(Mc>0);
[dif_pos(:,:,1),p_pos(:,:,1),val_pos(:,:,1),val_pos(:,:,2)] = SDL_Graph_Analyses(myfun,Mc(g11,:,:), Mc(g21,:,:)); % g1 vs g2
[dif_pos(:,:,2),p_pos(:,:,2),val_pos(:,:,3),val_pos(:,:,4)] = SDL_Graph_Analyses(myfun,Mc(g12,:,:), Mc(g22,:,:)); % g1 vs g3
[dif_pos(:,:,3),p_pos(:,:,3),val_pos(:,:,5),val_pos(:,:,6)] = SDL_Graph_Analyses(myfun,Mc(g13,:,:), Mc(g23,:,:)); % g1 vs g4
[dif_pos(:,:,4),p_pos(:,:,4),val_pos(:,:,7),val_pos(:,:,8)] = SDL_Graph_Analyses(myfun,Mc(g14,:,:), Mc(g24,:,:)); % g1 vs g5
%[dif_pos(:,:,5),p_pos(:,:,5),val_pos(:,:,9),val_pos(:,:,10)] = SDL_Graph_Analyses(myfun,Mc(g15,:,:), Mc(g25,:,:)); % g2 vs g3
%[dif_pos(:,:,6),p_pos(:,:,6),val_pos(:,:,11),val_pos(:,:,12)] = SDL_Graph_Analyses(myfun,Mc(g16,:,:), Mc(g26,:,:)); % g2 vs g4
%[dif_pos(:,:,7),p_pos(:,:,7),val_pos(:,:,13),val_pos(:,:,14)] = SDL_Graph_Analyses(myfun,Mc(g17,:,:), Mc(g27,:,:)); % g2 vs g5
[dif_pos(:,:,5),p_pos(:,:,5),val_pos(:,:,9),val_pos(:,:,10)] = SDL_Graph_Analyses(myfun,Mc(g18,:,:), Mc(g28,:,:)); % g3 vs g4

% neg values
Thr = Thrr(4*N+1:6*N,:); WC = WCC(4*N+1:6*N,:);
% ThC  = Thr(find(WC==max(WC(:,1)))); % the threshold corresponding to the maximum value of minimum wring cost
% Mc  = double(Mci<ThC(1));
Mc = SDL_threshold_proportional(-1*Mci,max(WC(:,1)));
Mc = double(Mc>0);
[dif_neg(:,:,1),p_neg(:,:,1),val_neg(:,:,1),val_neg(:,:,2)] = SDL_Graph_Analyses(myfun,Mc(g11,:,:), Mc(g21,:,:)); % g1 vs g2
[dif_neg(:,:,2),p_neg(:,:,2),val_neg(:,:,3),val_neg(:,:,4)] = SDL_Graph_Analyses(myfun,Mc(g12,:,:), Mc(g22,:,:)); % g1 vs g3
[dif_neg(:,:,3),p_neg(:,:,3),val_neg(:,:,5),val_neg(:,:,6)] = SDL_Graph_Analyses(myfun,Mc(g13,:,:), Mc(g23,:,:)); % g1 vs g4
[dif_neg(:,:,4),p_neg(:,:,4),val_neg(:,:,7),val_neg(:,:,8)] = SDL_Graph_Analyses(myfun,Mc(g14,:,:), Mc(g24,:,:)); % g1 vs g5
%[dif_neg(:,:,5),p_neg(:,:,5),val_neg(:,:,9),val_neg(:,:,10)] = SDL_Graph_Analyses(myfun,Mc(g15,:,:), Mc(g25,:,:)); % g2 vs g3
%[dif_neg(:,:,6),p_neg(:,:,6),val_neg(:,:,11),val_neg(:,:,12)] = SDL_Graph_Analyses(myfun,Mc(g16,:,:), Mc(g26,:,:)); % g2 vs g4
%[dif_neg(:,:,7),p_neg(:,:,7),val_neg(:,:,13),val_neg(:,:,14)] = SDL_Graph_Analyses(myfun,Mc(g17,:,:), Mc(g27,:,:)); % g2 vs g5
[dif_neg(:,:,5),p_neg(:,:,5),val_neg(:,:,9),val_neg(:,:,10)] = SDL_Graph_Analyses(myfun,Mc(g18,:,:), Mc(g28,:,:)); % g3 vs g4

end


% Graph theory comparisons between M1 and M2
% using Matlab handle function
function [dif, p, v1, v2] = SDL_Graph_Analyses(myfun,M1,M2)
% Input
% ------ myfun, Matlab handle function, e.g. @degrees_und
% ------ M1, binary matrix of group 1
% ------ M2, binary matrix of group 2
% Output
% ------ dif, between-group differences of nodal values
% ------ p,   p list of between-group differences across nodes
% ------ v1,  raw value of group1
% ------ v2,  raw value of group2

N = size(M1,3); % N = SDL.N + 1
if min(size(myfun(M1(:,:,1)))) == 1 % output is just 1 dimension
    for j = 1:N % per permutation
        D1(:,j) = myfun(M1(:,:,j)); % nodal centrality of raw data + randomized data, group1
        D2(:,j) = myfun(M2(:,:,j)); % nodal centrality of raw data + randomized data, group2
    end
    v1 = D1(:,1); % nodal centrality of raw data, group1
    v2 = D2(:,1); % nodal centrality of raw data, group2
    diff = D1 - D2; % M1 - M2 for all permutations, dim=2
else % output is more than 1 dimension
    for j = 1:N % per permutation
        D1(:,:,j) = myfun(M1(:,:,j));
        D2(:,:,j) = myfun(M2(:,:,j));
    end
    if strcmp(func2str(myfun),'distance_bin') % closeness = 1/distance
        DD1 = 1./D1; DD1(~isfinite(DD1))=0; %
        DD2 = 1./D2; DD2(~isfinite(DD2))=0; %
        diff = reshape(sum(DD1-DD2,2),size(DD1,1),size(DD1,3)); % M1 - M2 for all permutations, dim=2 after reshape
        v1 = sum(DD1(:,:,1),2);% nodal centrality of raw data, group1
        v2 = sum(DD2(:,:,1),2);% nodal centrality of raw data, group2
    else
        %         diff = reshape(sum(D1-D2,2),size(D1,1),size(D1,3)); % M1 - M2 for all permutations, dim=2 after reshape
        diff = D1 - D2; % M1 - M2 for all permutations, dim=3
        v1 = D1(:,:,1); % nodal centrality of raw data, group1
        v2 = D2(:,:,1); % nodal centrality of raw data, group2
    end
end

if ndims(diff) == 3 % dim=3
    dif  = diff(:,:,1);         % M1 - M2 for raw data
    p    = SDL_Pvalue(diff(:,:,1), diff(:,:,2:end), N-1); %p values, work for both dim 2 & 3
else % dim=2
    dif  = diff(:,1);         % M1 - M2 for raw data
    p    = SDL_Pvalue(diff(:,1), diff(:,2:end), N-1); %p values, work for both dim 2 & 3
end
end

% p value of between-group comparison
function P = SDL_Pvalue(Dreal0, DrandOrig0, M)
%==========================================================================
%   This function is based on gretna_permutation_test.m in GRETNA toolbox
% Inputs
% ------ Dreal0:     real value,
% ------ DrandOrig0: vector contains M times permutation values
% ------ M:         number of permutation
% Outputs
% --- P:      statistical p value

dimflag = 2; % dim of DrandOrig
[s1, s2] = size(Dreal0);
% reshape data if dim = 3
if length(size(DrandOrig0)) == 3 % dim = 3
    dimflag = 3;
    Dreal     = reshape(Dreal0,    numel(Dreal0),1);
    DrandOrig = reshape(DrandOrig0,numel(Dreal0),size(DrandOrig0,3));
else % dim = 2
    Dreal     = Dreal0;
    DrandOrig = DrandOrig0;
end

% p values calculation
for k = 1:numel(Dreal) % per node or edge
    
    Delta.real = Dreal(k);
    Delta.rand = DrandOrig(k,:);
    
    if Delta.real > 0
        P(k,1) = (1+length(find(Delta.rand >= Delta.real)))/(M+1);
    else
        P(k,1) = (1+length(find(Delta.rand <= Delta.real)))/(M+1);
    end
end

% reshape p values of dim = 3
if dimflag == 3 % dim = 3
    P = reshape(P, s1,s2);
end

end
