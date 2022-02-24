function SDL_ResultsReport_Med(SDL)
% make report tables & figures
% Input
% ------ SDL, all important parameters
% Output
% ------ SDL, with a newly added Report, a table report 
%                --- labels of nodes, 
%                --- graph values,
%                --- statistical values (unc.)
%                --- statistical values (corr.)

%% Load data
fraw = SDL.raw;
fout = SDL.out;
img = SDL.img; % whether to produce brain maps or not
fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
fn = fullfile(fdir,['Reports_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);

if strcmp(SDL.ana_type{1},'med') || strcmp(SDL.ana_type{1},'medcov') % mediation matrix
    ido = SDL.XYM{3}{1}; % index of areas serving as mediator (thus removed from tbl)
else
    ido = [];
end
fn
load(fn); fprintf('Loaded: Results report in %s\n',fn);
BNVpath = fullfile(SDL.BoxPath,'BrainNetViewer_20171031'); % add the path of BrainNetViewer

% get labels of all cortical areas
if isfield(SDL,'T')
    tbl = SDL.T.Properties.VariableNames(2:149); 
else
    Ta = readtable(fraw,'sheet',SDL.data_type{1}(1:2));
    tbl = Ta.Properties.VariableNames(2:149); 
end

%% Display results tables
%     T_Mc_abs  = SDL_Reports_Node(SDL.val_Mc_abs,SDL.dif_Mc_abs, SDL.p_Mc_abs, tbl,ido,fn,'Mc_abs',BNVpath,img) % abs correlations or regressions
%     T_Mc_pos  = SDL_Reports_Node(SDL.val_Mc_pos,SDL.dif_Mc_pos, SDL.p_Mc_pos, tbl,ido,fn,'Mc_pos',BNVpath,img) % pos correlations or regressions
%     T_Mc_neg  = SDL_Reports_Node(SDL.val_Mc_neg,SDL.dif_Mc_neg, SDL.p_Mc_neg, tbl,ido,fn,'Mc_neg',BNVpath,img) % neg correlations or regressions

id = 1; % Mc
val1(:,1,:) = SDL.Node_Degree.val1_pos(:,:,id); val1(:,2,:) = SDL.Node_Betweenness.val1_pos(:,:,id); val1(:,3,:) = SDL.Node_Closeness.val1_pos(:,:,id);
val2(:,1,:) = SDL.Node_Degree.val2_pos(:,:,id); val2(:,2,:) = SDL.Node_Betweenness.val2_pos(:,:,id); val2(:,3,:) = SDL.Node_Closeness.val2_pos(:,:,id);
dif(:,1,:) = SDL.Node_Degree.dif_pos(:,:,id); dif(:,2,:) = SDL.Node_Betweenness.dif_pos(:,:,id); dif(:,3,:) = SDL.Node_Closeness.dif_pos(:,:,id);
pli(:,1,:) = SDL.Node_Degree.p_pos(:,:,id );  pli(:,2,:) = SDL.Node_Betweenness.p_pos(:,:,id );  pli(:,3,:) = SDL.Node_Closeness.p_pos(:,:,id);
T_Mc_pos  = SDL_Reports_Node(val1,val2,dif, pli, tbl,ido,fn,'Mc_pos',BNVpath,img,'Mc') % pos correlations or regressions

id = 2; % Mc1
val1(:,1,:) = SDL.Node_Degree.val1_pos(:,:,id); val1(:,2,:) = SDL.Node_Betweenness.val1_pos(:,:,id); val1(:,3,:) = SDL.Node_Closeness.val1_pos(:,:,id);
val2(:,1,:) = SDL.Node_Degree.val2_pos(:,:,id); val2(:,2,:) = SDL.Node_Betweenness.val2_pos(:,:,id); val2(:,3,:) = SDL.Node_Closeness.val2_pos(:,:,id);
dif(:,1,:) = SDL.Node_Degree.dif_pos(:,:,id); dif(:,2,:) = SDL.Node_Betweenness.dif_pos(:,:,id); dif(:,3,:) = SDL.Node_Closeness.dif_pos(:,:,id);
pli(:,1,:) = SDL.Node_Degree.p_pos(:,:,id );  pli(:,2,:) = SDL.Node_Betweenness.p_pos(:,:,id );  pli(:,3,:) = SDL.Node_Closeness.p_pos(:,:,id);
T_Mc1_pos  = SDL_Reports_Node(val1,val2,dif, pli, tbl,ido,fn,'Mc1_pos',BNVpath,img,'Mc1') % pos correlations or regressions


id = 3; % Mab
val1(:,1,:) = SDL.Node_Degree.val1_pos(:,:,id); val1(:,2,:) = SDL.Node_Betweenness.val1_pos(:,:,id); val1(:,3,:) = SDL.Node_Closeness.val1_pos(:,:,id);
val2(:,1,:) = SDL.Node_Degree.val2_pos(:,:,id); val2(:,2,:) = SDL.Node_Betweenness.val2_pos(:,:,id); val2(:,3,:) = SDL.Node_Closeness.val2_pos(:,:,id);
dif(:,1,:) = SDL.Node_Degree.dif_pos(:,:,id); dif(:,2,:) = SDL.Node_Betweenness.dif_pos(:,:,id); dif(:,3,:) = SDL.Node_Closeness.dif_pos(:,:,id);
pli(:,1,:) = SDL.Node_Degree.p_pos(:,:,id );  pli(:,2,:) = SDL.Node_Betweenness.p_pos(:,:,id );  pli(:,3,:) = SDL.Node_Closeness.p_pos(:,:,id);
T_Mab_pos  = SDL_Reports_Node(val1,val2,dif, pli, tbl,ido,fn,'Mab_pos',BNVpath,img,'Mab') % pos correlations or regressions

%% Load data
fdir = fullfile(fout,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat'])
load(fn); fprintf('Loaded: Mediation Matrix in %s\n',fn);
[W.WC_abs, W.WC_pos, W.WC_neg]  % wiring costs

disp('OK')
% END
end

% % make results tables for edge measures
% function T = SDL_Reports_Edge(diff,plist,tbl,ido,fn,ftype,BNVpath,img) 
% % Input
% % ------ diff, list of between-group value differences across nodes
% % ------ plist, list of p values of between-group comparisons
% % ------ tbl, labels of all nodes
% % ------ ido, nodes not included because they are mediator
% % ------ fn, filename of the Results Report file
% % ------ ftype, 'Mc','Mc1','Mab' (the latter two for mediation analyses only) combined with 'abs' ,'pos' or 'neg'
% % ------ BNVpath, path of Brain Net Viewer toolbox
% % ------ img, produce brain maps if img = 'Brain Maps = YES'
% % Output
% % ------ T, table of significant outputs for Mc, Mc1 and Mab
% 
% T = []; 
% tbl1  = tbl; 
% if ~isempty(ido) % if there is a mediator
%     tbl1(ido) = []; % remove mediator from labels of nodes
% end
% [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(plist,0.05,'pdep','no'); % FDR corrected p values
% 
% [x,y] = find(plist<0.005); % too many connections when < 0.05, so default using , 0.005
% 
% 
% end


% make results tables for nodal measures
function T = SDL_Reports_Node(val1,val2, diff,plist,tbl,ido,fn,ftype,BNVpath,img,idname)
% Input
% ------ val1&2,  list of values per group
% ------ diff, list of between-group value differences across nodes
% ------ plist,list of p values of between-group comparisons
% ------ tbl, labels of all nodes
% ------ ido, nodes not included because they are mediator
% ------ fn, filename of the Results Report file
% ------ ftype, 'Mc','Mc1','Mab' (the latter two for mediation analyses only) combined with 'abs' ,'pos' or 'neg'
% ------ BNVpath, path of Brain Net Viewer toolbox
% ------ img, produce brain maps if img = 'Brain Maps = YES'
% ------ idname, name of between-group comparisons, e.g. 'g1-g2'
% Output
% ------ T, table of significant outputs for Mc, Mc1 and Mab
T = []; 
tbl1  = tbl; 
if ~isempty(ido) % if there is a mediator
    tbl1(ido) = []; % remove mediator from labels of nodes
end
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(plist,0.05,'pdep','no'); % FDR corrected p values

idx         = find(sum(plist<0.05,2)>2); % index of areas with at least three centralities are < 0.05
idx0        = idx; 
sz = length(tbl1)/2;
idx0(idx>sz)= idx(idx>sz)-sz; % label at each hemisphere

T_Mc  = array2table(idx0,'VariableNames',{'NumHem'}); T_Mc.Nodes   = tbl1(idx)'; % the suffix '_thickavg' does not necessarily mean CT

T_Mc.val_degree_g1   = val1(idx,1);
T_Mc.val_degree_g2   = val2(idx,1);
T_Mc.val_betweenness_g1   = val1(idx,2);
T_Mc.val_betweenness_g2   = val2(idx,2);
T_Mc.val_closeness_g1   = val1(idx,3);
T_Mc.val_closeness_g2   = val2(idx,3);


T_Mc.dif_degree      = diff(idx,1); % use values of degree centrality
T_Mc.dif_betweenness = diff(idx,2); % use values of betweenness centrality
T_Mc.dif_closeness   = diff(idx,3); % use values of closeness centrality

T_Mc.p_degree        = plist(idx,1);
T_Mc.p_betweenness   = plist(idx,2);
T_Mc.p_closeness     = plist(idx,3);

T_Mc.pFDR_min  = min(adj_p(idx,:),[],2); % min of all 4 p values


% display results
[FILEPATH,NAME,EXT] = fileparts(fn);% Brain Net Viewer Displaying (and saving) Areas Associated Significant Components
fout = fullfile(FILEPATH,'Fig',[NAME,idname,'_',ftype,'.png']); % image containing the nodes associated with signficant component
T = T_Mc;
fprintf('Above Reports displayed in\n %s\n',fout);

if strcmp(img,'Brain Maps = YES')
    BrainNet(BNVpath,idx,diff,ido,fout);
else
    fprintf('No brain maps are required to produce\n');
end

end

% Display the nodes associated with significant components through BrainNet Viewer
function BrainNet(BNVpath,idx_list,dif_list,ido,fout)
    [FILEPATH,NAME,EXT] = fileparts(fout);mkdir(FILEPATH);
    % make edge file for the output image
    fEdge = fullfile(FILEPATH,[NAME,'.edge']);
    dlmwrite(fEdge,zeros(148,148),'delimiter',' ');
    % load sample node file
    fid = fopen('Node_sample.node');
    fline = textscan(fid,'%f%f%f%d%d%s\n'); % col 1-3, coordinates; 4, color; 5, size; 6, name
    fclose(fid);
    % revising node information
    fline{4} = zeros(148,1); fline{5} = zeros(148,1);
    
    if ~isempty(ido) % if there is a mediator
        for i = 1:6 % per column in fline
            for j = ido % per area as mediator
                fline{i}(j) = [];
            end
        end
    end
    
    if ~isempty(idx_list) % if there is any significant node
        for ic = 1:size(idx_list,1) % per significant node
            if dif_list(idx_list(ic),1) > 0 % g1 > g2, based on degree centrality, but should check whether consistent with other centralities.
                fline{4}(idx_list(ic)) = 2; % red color
            else % g1 < g2
                fline{4}(idx_list(ic)) = 0; % blue color
            end
            fline{5}(idx_list(ic)) = 1;  % change size
        end
    end

    % make node file for the output image
    fNode = fullfile(FILEPATH,[NAME,'.node']);
    fid = fopen(fNode,'w');
    for i=1:length(fline{1})
        fprintf(fid,'%1.2f\t%1.2f\t%1.2f\t%d\t%d\t%s\n',...
            fline{1}(i),fline{2}(i),fline{3}(i),fline{4}(i),fline{5}(i),fline{6}{i});
    end
    fclose(fid);
    
    % BrainNet Viewe Display
    warning('OFF'); rmpath(genpath(BNVpath)); warning('ON'); addpath(genpath(BNVpath)); % remove and then add path
    fnv   = fullfile(BNVpath,'Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv');
    fCfg  = 'Cfg.mat';
    fout  = fullfile(FILEPATH,[NAME,'.png']);
%     figure;
%     BrainNet_MapCfg(fnv,fNode,fEdge, fout,fCfg);
    BrainNet_MapCfg('/Users/ds366/Documents/MATLAB/SDLToolBox/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
    fNode,fEdge,fout,'/Volumes/dusom_morey/Data/Lab/Delin/Projects/ENIGMA_SCA_Rakesh/Scripts/Cfg.mat');
end
