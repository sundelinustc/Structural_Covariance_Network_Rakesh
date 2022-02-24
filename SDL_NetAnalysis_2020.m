                              %==============================================
% This script is for structural covariance analyses on cortical thickness or surface area networks,
% by Dr. Delin Sun, in Morey Lab, Duke University, 01-22-2020
%==============================================
clear;
close all;   

%% Parameters11
[FILEPATH,NAME,EXT] = fileparts(pwd);
SDL.path = FILEPATH;
% SDL.path = 'Z:/Data/Lab/Delin/Projects/ENIGMA_SCA_Rakesh';
% SDL.path = '/Volumes/dusom_morey-1/Data/Lab/Delin/Projects/ENIGMA_SCA_Rakesh'; % data stored at my Mac desktop, remotly access using X32
SDL.BoxPath = fullfile(SDL.path,'Scripts','BCT'); % path for most of the toolboxes
SDL.raw = fullfile(SDL.path,'Original','PGC_Data_Covariates_check_for_delin_withCAPS_original.xlsx'); % raw data file
SDL.out = fullfile(SDL.path,'Outputs');
   
% data type & analyses type
% CT = cortical thickness, SA = surface area
% corr = correlation, partialcorr = partial correlation,
% regr = regression, regrcov = regression with covariates (time consuming)
% med = mediation, medcov = mediation with covariates (time consuming)
% und = undirected connection, dir = directed connection
% ---data_type, ana_type, connection with direction or not, predefined X, Y and M
Ana = {
         'CT',           'corr', 'und', {},{},{},'';
         'SA',           'corr', 'und', {},{},{},'';
%     
%          'CT_PTSDsev10', 'corr', 'und', {},{},{},'';
%          'SA_PTSDsev10', 'corr', 'und', {},{},{},'';
%     
%          'CT_Dep',       'corr', 'und', {},{},{},'';
%          'SA_Dep',       'corr', 'und', {},{},{},'';
%     
%          'CT_Gender',    'corr', 'und', {},{},{},'';
%          'SA_Gender',    'corr', 'und', {},{},{},'';
    
%          'CT_Age10',     'corr', 'und', {},{},{},'';
%          'SA_Age10',     'corr', 'und', {},{},{},'';
    
    
    
    %     'CT', 'partialcorr', 'und', {},{},{},'';
    %     'SA', 'partialcorr', 'und', {},{},{},'';
    %
    %     'CT', 'regr',        'und', {},{},{},'';
    %     'SA', 'regr',        'und', {},{},{},'';
    %
    %'CT', 'med',        'dir', {},{},{[[6,7],74+[6,7]]},'X=[],Y=[],M=ACC'; % M = SN-ACC
    %'SA', 'med',        'dir', {},{},{[[6,7],74+[6,7]]},'X=[],Y=[],M=ACC'; % M = SN-ACC
    %'CT', 'med',        'dir', {},{},{[18,40,47,49]},   'X=[],Y=[],M=L aIns'; % M = SN-L SN
    %'SA', 'med',        'dir', {},{},{[18,40,47,49]},   'X=[],Y=[],M=L aIns'; % M = SN-L SN
    %'CT', 'med',        'dir', {},{},{74+[18,40,47,49]},'X=[],Y=[],M=R aIns'; % M = SN-R SN
    %'SA', 'med',        'dir', {},{},{74+[18,40,47,49]},'X=[],Y=[],M=R aIns'; % M = SN-R SN
    
    %     'CT', 'med',        'dir', {},{},{31,32,63,70,74+31,74+32,74+63,74+70},'X=[],Y=[],M=vmPFC'; % M = DMN-vmPFC
    %     'SA', 'med',        'dir', {},{},{31,32,63,70,74+31,74+32,74+63,74+70},'X=[],Y=[],M=vmPFC'; % M = DMN-vmPFC
    %     'CT', 'med',        'dir', {},{},{9,10,74+9,74+10},'X=[],Y=[],M=PCC'; % M = DMN-PCC
    %     'SA', 'med',        'dir', {},{},{9,10,74+9,74+10},'X=[],Y=[],M=PCC'; % M = DMN-PCC
    %    'CT', 'med',        'dir', {},{},{25,26,36,41,55},'X=[],Y=[],M=L TPJ'; % M = DMN-L TPJ
    %     'SA', 'med',        'dir', {},{},{25,26,36,41,55},'X=[],Y=[],M=L TPJ'; % M = DMN-L TPJ
    %    'CT', 'med',        'dir', {},{},{74+25,74+26,74+36,74+41,74+55},'X=[],Y=[],M=R TPJ'; % M = DMN-R TPJ
    %     'SA', 'med',        'dir', {},{},{74+25,74+26,74+36,74+41,74+55},'X=[],Y=[],M=R TPJ'; % M = DMN-R TPJ
    
    %     'CT', 'med',        'dir', {},{},{15,29,52,68},'X=[],Y=[],M=L DLPFC'; % M = FP-L DLPFC
    %     'SA', 'med',        'dir', {},{},{15,29,52,68},'X=[],Y=[],M=L DLPFC'; % M = FP-L DLPFC
    %     'CT', 'med',        'dir', {},{},{74+15,74+29,74+52,74+68},'X=[],Y=[],M=R DLPFC'; % M = FP-R DLPFC
    %     'SA', 'med',        'dir', {},{},{74+15,74+29,74+52,74+68},'X=[],Y=[],M=R DLPFC'; % M = FP-R DLPFC
    %     'CT', 'med',        'dir', {},{},{27,28,56,67},'X=[],Y=[],M=L PP'; % M = FP-L PP
    %     'SA', 'med',        'dir', {},{},{27,28,56,67},'X=[],Y=[],M=L PP'; % M = FP-L PP
    %     'CT', 'med',        'dir', {},{},{74+27,74+28,74+56,74+67},'X=[],Y=[],M=R PP'; % M = FP-R PP
    %     'SA', 'med',        'dir', {},{},{74+27,74+28,74+56,74+67},'X=[],Y=[],M=R PP'; % M = FP-R PP
    
    };

SDL.group = {'PTSD','CONT'};
SDL.N = 5000; % number of permutation, 10,000 for publication purpose

disp (SDL);


% %% Analyses
for i = 1:size(Ana,1) % per pre-defined analysis
    % Types for analyses
    SDL.data_type = Ana(i,1);
    SDL.ana_type  = Ana(i,2);
    SDL.XYM       = Ana(i,4:7); % for mediation analyses only
    %SDL.img       = 'Brain Maps = YES'; % make figures of brain maps about between-group comparison results
    SDL.img       = 'Brain Maps = NO'; % make figures of brain maps about between-group comparison results
    
    
%     %load and clean raw data
%      SDL_Load_Clean(SDL); % load and clean raw data
%      SDL_ComBat(SDL); % Harmonization data across study sites using ComBat
%      SDL_LME_residules(SDL); % residules after LME to exclude age, age^2, sex, and mean CT or SA
    
    if strcmp(SDL.ana_type{1},'corr') || strcmp(SDL.ana_type{1},'corr') % correlation matrix
        % Network Generation and Graph Analyses
        if strcmp(SDL.data_type{1},'CT') || strcmp(SDL.data_type{1},'SA') % main effects of PTSD diagnosis
%                     SDL_matrix(SDL); % calculate correlation, regression or mediation matrix (C,C1,ab) per permutation
%                    SDL_WiringCost_Th(SDL); % minimum wiring cost, and the mediation matrix after thresholding
                    SDL_Net_Analyses_NoMed(SDL); % Graph theory analyses, non-mediation matrix
%            SDL_ResultsReport(SDL); % report findings
       %elseif strcmp(SDL.data_type{1},'CT_PTSDsev10') || strcmp(SDL.data_type{1},'SA_PTSDsev10') % effects of PTSD symptom syverity
        %             SDL_matrix_PTSDsev10(SDL); % calculate correlation, regression or mediation matrix (C,C1,ab) per permutation
         %            SDL_WiringCost_Th_PTSDsev10(SDL); % minimum wiring cost, and the mediation matrix after thresholding
          %           SDL_Net_Analyses_NoMed_PTSDsev10(SDL); % Graph theory analyses, non-mediation matrix
           % SDL_ResultsReport_PTSDsev10(SDL); % report findings
%             if strcmp(SDL.data_type{1},'CT_Age10') || strcmp(SDL.data_type{1},'SA_Age10')
%                      SDL_matrix_Age10(SDL); % calculate correlation, regression or mediation matrix (C,C1,ab) per permutation
%                      SDL_WiringCost_Th_Age10(SDL); % minimum wiring cost, and the mediation matrix after thresholding
%                      SDL_Net_Analyses_NoMed_Age10(SDL); % Graph theory analyses, non-mediation matrix
%             SDL_ResultsReport_Age10(SDL); % report findings
            
            %
      %  else including PTSD x Gender and PTSD x Depression
                    % SDL_matrix_All4(SDL); % calculate correlation, regression or mediation matrix (C,C1,ab) per permutation
                     %SDL_WiringCost_Th_All4(SDL); % minimum wiring cost, and the mediation matrix after thresholding
                     %SDL_Net_Analyses_NoMed_All4(SDL); % Graph theory analyses, non-mediation matrix
                     %SDL_ResultsReport_All4(SDL); % report findings
            
        end
   % else  mediatin analyses
    %     SDL_matrix(SDL); % calculate correlation, regression or mediation matrix (C,C1,ab) per permutation
     %     SDL_WiringCost_Th(SDL); % minimum wiring cost, and the mediation matrix after thresholding
      %   SDL_Net_Analyses(SDL);
       % SDL_ResultsReport_Med(SDL); % report findings
    end
    
    %     SDL_extra(SDL); % for some extra analyses
end
