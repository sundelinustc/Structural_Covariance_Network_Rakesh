function SDL_extra(SDL)

% for some extra analyses





%% (1) Between-group difference of CT/SA
% load data
fn = fullfile(SDL.path,'Outputs','CleanData');
if ~strcmp(SDL.ana_type{1},'med') % correlation or regression matrix
    fn = fullfile(fn,['Data_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
else
    fn = fullfile(fn,['Data_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
end
load(fn);

% Analyses
T = SDL.T; Coef = []; plist = [];
roi_list = cell2mat(SDL.XYM{3});
for j = roi_list % per column of data values of interest
    tbl = T(:,{'Group','Age','Gender','Site'}); tbl.data = T{:,j};
    lme = fitlme(tbl,'data ~ Group + Age + Gender + (1|Site)');
    Coef  = [Coef,lme.Coefficients{2,2}];  % lme regression coefficient for "Group_PTSD"
    plist = [plist,lme.Coefficients{2,6}]; % lme regression p value     for "Group_PTSD"
end
T_show = array2table(roi_list','VariableNames',{'ROI_No'});
T_show.ROI_Name = SDL.T.Properties.VariableNames(1+roi_list)';
T_show.Coef = Coef';
T_show.p_value = plist';
T_show

disp('OK');





%% End
end