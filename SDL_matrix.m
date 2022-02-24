function SDL_matrix(SDL)

%==============================================
% (1) remove effects of covariates
% (2) permutation N-1 times by shuffling group labels
% (3) calculate mediation matrx for each permutation
% Input
% ------ SDL, variable contains all the important path information
% Output
% ------ SDL, 
%        MMc,   total effect
%        MMc1,  direct effect
%        MMab,  indirect effect
%        They are also saved in 'Results_MediationMatrix_xxxx.mat' in Results folder 
%==============================================


fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn,'T');
R = T{:,2:149}; % residules

% %% (1) remove effects of covariates by using residuals of lme model
% T = SDL.T; % the input raw dataset after cleaning (each row-subject; each
% %           column-measure), n=size(T,1)
% if strcmp(SDL.ana_type{1},'med') % mediaton matrix
%     % residules of mediator
%     mlist = SDL.XYM{3}{:}; % index of areas serve as mediator
%     tbl = T(:,{'Age','Gender','Site'}); tbl.data = mean(T{:,mlist},2); % mean values of areas as mediator
%     lme = fitlme(tbl,'data ~ Age + Gender + (1|Site)');
%     RMV = residuals(lme); % residules of mediator
% else
%     mlist = [];
% end
% 
% % residules of other areas   
% R = []; % matrix containing residuals after lme
% k = 0;
% for j = 2:149 % per column of data values
%     if any(mlist==(j-1))
%     else
%         k = k + 1;
%         tbl = T(:,{'Age','Gender','Site'}); tbl.data = T{:,j};
%         lme = fitlme(tbl,'data ~ Age + Gender + (1|Site)');
%         R(:,k) = residuals(lme);
%     end
% end
% fprintf('Completed: Remove effects of covariates by using residuals of lme model\n');


%% (2) permutation N times by shuffling group labels
g1 = strcmp(T.Group,SDL.group{1}); g2 = strcmp(T.Group,SDL.group{2}); % index of group 1 & 2

N = SDL.N; % times of permutation: N>1, calculate the raw mediation matrix,
%           as well as N-1 mediation matrix with shuffling labels between
%           groups

if strcmp(SDL.ana_type{1},'corr') || strcmp(SDL.ana_type{1},'partialcorr') % correlation/partial-correlation matrix
    if strcmp(SDL.ana_type{1},'corr')
        SDL_corr = @corr;        % linear partial correlation coefficients between pairs of variables without controlling other variables
    else
        SDL_corr = @partialcorr; % linear partial correlation coefficients between pairs of variables with controlling other variables
    end
    for i = 1:N % per permutation
        tic;
        if i == 1 % no permutation
            RHO1 = SDL_corr(R(g1,:)); % linear partial correlation coefficients between pairs of variables in group 1
            RHO2 = SDL_corr(R(g2,:)); % linear partial correlation coefficients between pairs of variables in group 2
        else % permutation by shuffling group labels
            rmlist = randperm(sum(g1)+sum(g2));
            RHO1 = SDL_corr(R(rmlist(g1),:)); % linear partial correlation coefficients between pairs of variables in group 1
            RHO2 = SDL_corr(R(rmlist(g2),:)); % linear partial correlation coefficients between pairs of variables in group 2
        end
        MRHO(:,:,i) = [RHO1;RHO2];
        fprintf('Completed: permutation %d/%d for correlation matrix\t',i,N); toc;
    end
    fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
    fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
%     SDL.CorrMatrix = fn; % store the name of Results_MediationMatrix_xxx.mat file
    save(fn,'-v7.3'); fprintf('Saved: Matrix in %s\n',fdir);
    
elseif strcmp(SDL.ana_type{1},'regr') % regression matrix
    for i = 1:N % per permutation
        tic;
        if i == 1 % no permutation
            B1 = SDL_regression(R(g1,:)); % regression coefficients between pairs of variables in group 1
            B2 = SDL_regression(R(g2,:)); % regression coefficients between pairs of variables in group 1
        else % permutation by shuffling group labels
            rmlist = randperm(sum(g1)+sum(g2));
            B1 = SDL_regression(R(rmlist(g1),:)); % regression coefficients between pairs of variables in group 1
            B2 = SDL_regression(R(rmlist(g2),:)); % regression coefficients between pairs of variables in group 1
        end
%         MB(:,:,i) = [B1;B2]; 
        MB(:,:,i) = [(B1+B1')/2;(B2+B2')/2]; % to make the matrix symmetric
        fprintf('Completed: permutation %d/%d\t',i,N); toc;
    end
    fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
    fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
    SDL.RegrMatrix = fn; % store the name of Results_MediationMatrix_xxx.mat file
    save(fn,'-v7.3'); fprintf('Saved: Regression Matrix in %s\n',fn);
    
elseif strcmp(SDL.ana_type{1},'med') % mediation matrix
    Mlist = SDL.XYM{3}{:}; % the region number of mediators
    RMV = R(:,Mlist); RMV = mean(RMV,2); % mean value of mediators
    R(:,Mlist) = []; % data after removing the mediators
    for i = 1:N % per permutation
        tic;
        if i == 1 % no permutation
            [M1c,M1c1,M1ab] = SDL_mdiation(R(g1,:),RMV(g1,:)); % mediation matrix for group 1
            [M2c,M2c1,M2ab] = SDL_mdiation(R(g2,:),RMV(g2,:)); % mediation matrix for group 2
        else % permutation by shuffling group labels
            rmlist = randperm(sum(g1)+sum(g2));
            [M1c,M1c1,M1ab] = SDL_mdiation(R(rmlist(g1),:),RMV(rmlist(g1),:)); % mediation matrix for group 1
            [M2c,M2c1,M2ab] = SDL_mdiation(R(rmlist(g2),:),RMV(rmlist(g2),:)); % mediation matrix for group 2
        end
        MMc(:,:,i) = [M1c;M2c]; MMc1(:,:,i) = [M1c1;M2c1]; MMab(:,:,i) = [M1ab;M2ab];
        fprintf('Completed: permutation %d/%d\t',i,N); toc;
    end
    fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
    fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
    SDL.MedMatrix = fn; % store the name of Results_MediationMatrix_xxx.mat file
    save(fn,'-v7.3'); fprintf('Saved: Mediation Matrix in %s\n',fn);
    
else
end
disp('Matrix OK!');
% end
end

% calculate the regression matrix
function B = SDL_regression(D)
% Input
% ------ R, matrix of residuals
% Output
% ------ B, regression matrix

for i = 1:size(D,2)
    for j = 1:size(D,2)
        if i == j % diagnal = 0
            B(i,j) = 0;
        else
            X = D(:,i); Y = D(:,j); 
%             Cov = D; Cov(:,[i,j]) = []; % Covariates are too time consuming (~200 sec per permutation)
            BB = regress(Y,[ones(size(X)),X]); 
            B(i,j) = BB(2); % 1st coefficient is for intercept, 2nd is for X
        end
    end
end

end


% calculate the mediation matrix
function [Mc,Mc1,Mab] = SDL_mdiation(R,RMV)
% Input
% ------ R, matrix of residuals
% ------ RMV, residules of mediator
% Output
% ------ M, mediation matrix

M = RMV; D = R;
% Calculate C, C1 and a*b using regression analyses (instead of mediation toolbox, which is too slow)
for i = 1:size(D,2) % per X region
    for j = 1:size(D,2) % per Y region
        if i == j % diagnal = 0
            Mc(i,j) = 0; Mc1(i,j) = 0; Mab(i,j) = 0;
        else
            X = D(:,i); Y = D(:,j); 
%             Cov = D; Cov(:,[i,j]) = []; % Covariates are too time consuming (~200 sec per permutation)
            B = regress(Y,[ones(size(X)),X]); 
            Mc(i,j)  = B(2); % all relationship between X and Y
            B = regress(Y,[ones(size(X)),X,M]); 
            Mc1(i,j) = B(2); % direct (unmediated) relationship between X and Y
            Mab(i,j) = Mc(i,j) - Mc1(i,j); % mediated relationship between X and Y
        end
    end
end
% % Below, Mc1 and Mab are not raw values but proportion of Mc
% Mc1 = Mc1./Mc; Mc1 = Mc1/max(max(abs(Mc1))); % normalize values between -1 and 1
% Mab = Mab./Mc; Mab = Mab/max(max(abs(Mab))); % normalize values between -1 and 1
end