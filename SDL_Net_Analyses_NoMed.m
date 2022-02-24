  function SDL_Net_Analyses_NoMed(SDL)

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
fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Matrix in %s\n',fn);
Thr = [W.Thr_abs;W.Thr_pos;W.Thr_neg]; % Thresholds
WC  = [W.WC_abs; W.WC_pos; W.WC_neg];  % wiring costs


if strcmp(SDL.ana_type{1},'corr') ||strcmp(SDL.ana_type{1},'partialcorr') % correlation matrix
    MMc = MRHO;
elseif strcmp(SDL.ana_type{1},'regr') || strcmp(SDL.ana_type{1},'regrcov') % regression matrix
    MMc = MB;
else
end
atype = 'und'; % connection without direction

%% Nodal Measures
% (1) Nodal Degree Centrality
% (1) Nodal Degree Centrality
% [dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@degrees_und,MMc,WC,Thr);
% SDL.Node_Degree =      struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);
% % (2) Nodal Betweenness Centrality
% [dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@betweenness_bin,MMc,WC,Thr);
% SDL.Node_Betweenness = struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);
% % (3) Nodal Closeness Centrality (i.e. the reciprocal of distance, i.e. close_g1vg2 = 1/distance1-1/distance2)
% [dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@distance_bin,MMc,WC,Thr);
% SDL.Node_Closeness =   struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);
% % (4) Eigenvector Centrality
% [dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@eigenvector_centrality_und,MMc,WC,Thr);
% SDL.Node_Eigenvector = struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);
% (5) Community_louvain
[dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@community_louvain,MMc,WC,Thr);
SDL.community_louvain =struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);
[dif_abs,p_abs,val_abs,dif_pos,p_pos,val_pos,dif_neg,p_neg,val_neg] = SDL_Graph_output(@modularity_und,MMc,WC,Thr);
SDL.modularity       = struct('dif_abs',dif_abs,'p_abs',p_abs,'val_abs',val_abs,'dif_pos',dif_pos,'p_pos',p_pos,'val_pos',val_pos,'dif_neg',dif_neg,'p_neg',p_neg,'val_neg',val_neg);


% % Mc_abs
% SDL.val_Mc_abs = [SDL.Node_Degree.val_abs(:,1,:),SDL.Node_Betweenness.val_abs(:,1,:),SDL.Node_Closeness.val_abs(:,1,:),SDL.Node_Eigenvector.val_abs(:,1,:)]; 
% SDL.dif_Mc_abs = [SDL.Node_Degree.dif_abs(:,1),SDL.Node_Betweenness.dif_abs(:,1),SDL.Node_Closeness.dif_abs(:,1),SDL.Node_Eigenvector.dif_abs(:,1)]; 
% SDL.p_Mc_abs   = [SDL.Node_Degree.p_abs(:,1),  SDL.Node_Betweenness.p_abs(:,1),  SDL.Node_Closeness.p_abs(:,1),  SDL.Node_Eigenvector.p_abs(:,1)];
% 
% % Mc_pos
% SDL.val_Mc_pos = [SDL.Node_Degree.val_pos(:,1,:),SDL.Node_Betweenness.val_pos(:,1,:),SDL.Node_Closeness.val_pos(:,1,:),SDL.Node_Eigenvector.val_pos(:,1,:)]; 
% SDL.dif_Mc_pos = [SDL.Node_Degree.dif_pos(:,1),SDL.Node_Betweenness.dif_pos(:,1),SDL.Node_Closeness.dif_pos(:,1),SDL.Node_Eigenvector.dif_pos(:,1)]; 
% SDL.p_Mc_pos   = [SDL.Node_Degree.p_pos(:,1),  SDL.Node_Betweenness.p_pos(:,1),  SDL.Node_Closeness.p_pos(:,1),  SDL.Node_Eigenvector.p_pos(:,1)];
% 
% % Mc_neg
% SDL.val_Mc_neg = [SDL.Node_Degree.val_neg(:,1,:),SDL.Node_Betweenness.val_neg(:,1,:),SDL.Node_Closeness.val_neg(:,1,:),SDL.Node_Eigenvector.val_neg(:,1,:)]; 
% SDL.dif_Mc_neg = [SDL.Node_Degree.dif_neg(:,1),SDL.Node_Betweenness.dif_neg(:,1),SDL.Node_Closeness.dif_neg(:,1),SDL.Node_Eigenvector.dif_neg(:,1)]; 
% SDL.p_Mc_neg   = [SDL.Node_Degree.p_neg(:,1),  SDL.Node_Betweenness.p_neg(:,1),  SDL.Node_Closeness.p_neg(:,1),  SDL.Node_Eigenvector.p_neg(:,1)];

% %% Edge Measures
% % (1) Edge Betweenness Centrality
% [dif_abs,p_abs,dif_pos,p_pos,dif_neg,p_neg] = SDL_Graph_output(@edge_betweenness_bin,MMc,WC,Thr);
% SDL.Edge_Betweenness = struct('dif_abs',dif_abs,'p_abs',p_abs,'dif_pos',dif_pos,'p_pos',p_pos,'dif_neg',dif_neg,'p_neg',p_neg);

%% save results
fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
fn = fullfile(fdir,['Reports_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'com_mod.mat']);
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

k = 0.5 * size(Mci,1); g1 = [1:k]; g2 = [k+1:2*k]; % index of two groups

% absolute values
Thr = Thrr(1:2,:); WC = WCC(1:2,:);
% ThC  = Thr(find(WC==max(WC(:,1)))); % the threshold corresponding to the maximum value of minimum wring cost
Mc = SDL_threshold_proportional(abs(Mci),max(WC(:,1)));
Mc = double(Mc>0);
[dif_abs(:,:,1),p_abs(:,:,1),val_abs(:,:,1),val_abs(:,:,2)] = SDL_Graph_Analyses(myfun,Mc(g1,:,:), Mc(g2,:,:));

% pos values
Thr = Thrr(3:4,:); WC = WCC(3:4,:);
% ThC  = Thr(find(WC==max(WC(:,1)))); % the threshold corresponding to the maximum value of minimum wring cost
% Mc  = double(Mci>ThC(1));
Mc = SDL_threshold_proportional(Mci,max(WC(:,1)));
Mc = double(Mc>0);
[dif_pos(:,:,1),p_pos(:,:,1),val_pos(:,:,1),val_pos(:,:,2)] = SDL_Graph_Analyses(myfun,Mc(g1,:,:), Mc(g2,:,:));

% neg values
Thr = Thrr(5:6,:); WC = WCC(5:6,:);
% ThC  = Thr(find(WC==max(WC(:,1)))); % the threshold corresponding to the maximum value of minimum wring cost
% Mc  = double(Mci<ThC(1));
Mc = SDL_threshold_proportional(-1*Mci,max(WC(:,1)));
Mc = double(Mc>0);
[dif_neg(:,:,1),p_neg(:,:,1),val_neg(:,:,1),val_neg(:,:,2)] = SDL_Graph_Analyses(myfun,Mc(g1,:,:), Mc(g2,:,:));

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
% ------ p, p list of between-group differences across nodes
% ------ v1,  raw value of group1
% ------ v2,  raw value of group2
Nodes_of_interest = [3,29,33,34,35,36,45,51,55,68,69,73,77,103,107,108,109,110,119,125,129,142,143,147];
M1 = M1(Nodes_of_interest,Nodes_of_interest,:);
M2 = M2(Nodes_of_interest,Nodes_of_interest,:); 

N = size(M1,3); % N = SDL.N + 1
if strcmp(func2str(myfun),'community_louvain') || strcmp(func2str(myfun),'modularity_und')
    k = 1;
    for j = 1:N % per permutation
        try
            [~, D1(:,k)] = myfun(M1(:,:,j));
            [~, D2(:,k)] = myfun(M2(:,:,j));
            k = k + 1;
        catch
        end
    end
    diff = D1 - D2; % M1 - M2 for all permutations, dim=2
    v1 = D1(:,:,1); % nodal centrality of raw data, group1
    v2 = D2(:,:,1); % nodal centrality of raw data, group2
else
    if min(size(myfun(M1(:,:,1)))) == 1 % output is just 1 dimension
        for j = 1:N % per permutation
            D1(:,j) = myfun(M1(:,:,j));
            D2(:,j) = myfun(M2(:,:,j));
        end
        diff = D1 - D2; % M1 - M2 for all permutations, dim=2
        v1 = D1(:,:,1); % nodal centrality of raw data, group1
        v2 = D2(:,:,1); % nodal centrality of raw data, group2
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
            % diff = reshape(sum(D1-D2,2),size(D1,1),size(D1,3)); % M1 - M2 for all permutations, dim=2 after reshape
            diff = D1 - D2; % M1 - M2 for all permutations, dim=3
            v1 = D1(:,:,1); % nodal centrality of raw data, group1
            v2 = D2(:,:,1); % nodal centrality of raw data, group2
        end
    end
end

if ndims(diff) == 3 % dim=3
    N = size(diff,3);
    dif  = diff(:,:,1);         % M1 - M2 for raw data
    p    = SDL_Pvalue(diff(:,:,1), diff(:,:,2:end), N-1); %p values, work for both dim 2 & 3
else % dim=2
    N = size(diff,2);
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
