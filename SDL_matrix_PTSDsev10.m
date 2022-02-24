function SDL_matrix_PTSDsev10(SDL)

%==============================================
% (1) permutation N-1 times by shuffling group labels
% (2) calculate correlation matrx for each permutation
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
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
T = T(T.CAPStype==4,:); % only subjects with CAPS-IV are included
R = T{:,2:149};

a = prctile(T.CAPS,[20 40 60 80]) % divide into 5 groups with compariable sample size
a = [4 20 46 68] % based on SA data, CT data are similiar
[sum(T.CAPS<a(1)),...
    sum(T.CAPS>=a(1)&T.CAPS<a(2)),...
    sum(T.CAPS>=a(2)&T.CAPS<a(3)),...
    sum(T.CAPS>=a(3)&T.CAPS<a(4)),...
    sum(T.CAPS>=a(4))]

idx1 = T.CAPS< a(1)               & T.CAPStype==4; [idx1,~] = find(idx1==1);   % group 1
idx2 = T.CAPS>=a(1) & T.CAPS<a(2) & T.CAPStype==4; [idx2,~] = find(idx2==1);   % group 2
idx3 = T.CAPS>=a(2) & T.CAPS<a(3) & T.CAPStype==4; [idx3,~] = find(idx3==1);   % group 3
idx4 = T.CAPS>=a(3) & T.CAPS<a(4) & T.CAPStype==4; [idx4,~] = find(idx4==1);   % group 4
idx5 = T.CAPS>=a(4) & T.CAPS<=120 & T.CAPStype==4; [idx5,~] = find(idx5==1);   % group 5

N1 = length(idx1);
N2 = length([idx1;idx2]);
N3 = length([idx1;idx2;idx3]);
N4 = length([idx1;idx2;idx3;idx4]);
N5 = length([idx1;idx2;idx3;idx4;idx5]);


%% (2) permutation N times by shuffling group labels
% g1 = strcmp(T.Group,SDL.group{1}); g2 = strcmp(T.Group,SDL.group{2}); % index of group 1 & 2

N = SDL.N; % times of permutation: N=0, calculate the raw mediation matrix;
%           N>1, N times correlation, partial correlation, or mediation matrix with shuffling labels between groups

if strcmp(SDL.ana_type{1},'corr') || strcmp(SDL.ana_type{1},'partialcorr') % correlation/partial-correlation matrix
    if strcmp(SDL.ana_type{1},'corr')
        SDL_corr = @corr;        % linear partial correlation coefficients between pairs of variables without controlling other variables
    else
        SDL_corr = @partialcorr; % linear partial correlation coefficients between pairs of variables with controlling other variables
    end
    for i = 1:N+1 % per permutation
        tic;
        if i == 1 % no permutation
            RHO1 = SDL_corr(R(idx1,:)); % linear (partial) correlation coefficients between pairs of variables in group 1
            RHO2 = SDL_corr(R(idx2,:)); % linear (partial) correlation coefficients between pairs of variables in group 2
            RHO3 = SDL_corr(R(idx3,:)); % linear (partial) correlation coefficients between pairs of variables in group 3
            RHO4 = SDL_corr(R(idx4,:)); % linear (partial) correlation coefficients between pairs of variables in group 4
            RHO5 = SDL_corr(R(idx5,:)); % linear (partial) correlation coefficients between pairs of variables in group 5
        else % permutation by shuffling group labels
            rmlist = [idx1;idx2;idx3;idx4;idx5]; % combining all subjects
            rl = randperm(length(rmlist)); % shuffle them
            RHO1 = SDL_corr(R(rmlist(rl(   1:N1)),:)); % linear partial correlation coefficients between pairs of variables in group 1
            RHO2 = SDL_corr(R(rmlist(rl(N1+1:N2)),:)); % linear partial correlation coefficients between pairs of variables in group 2
            RHO3 = SDL_corr(R(rmlist(rl(N2+1:N3)),:)); % linear partial correlation coefficients between pairs of variables in group 3
            RHO4 = SDL_corr(R(rmlist(rl(N3+1:N4)),:)); % linear partial correlation coefficients between pairs of variables in group 4
            RHO5 = SDL_corr(R(rmlist(rl(N4+1:N5)),:)); % linear partial correlation coefficients between pairs of variables in group 5
        end
        MRHO(:,:,i) = [RHO1;RHO2;RHO3;RHO4;RHO5];
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
for i = 1:size(D,2)
    for j = 1:size(D,2)
        if i == j % diagnal = 0
            Mc(i,j) = 0; Mc1(i,j) = 0; Mab(i,j) = 0;
        else
            X = D(:,i); Y = D(:,j); 
%             Cov = D; Cov(:,[i,j]) = []; % Covariates are too time consuming (~200 sec per permutation)
            B = regress(Y,[ones(size(X)),X]); 
            Mc(i,j)  = B(2);
            B = regress(Y,[ones(size(X)),X,M]); 
            Mc1(i,j) = B(2); 
            Mab(i,j) = Mc(i,j) - Mc1(i,j);
        end
    end
end
% % Below, Mc1 and Mab are not raw values but proportion of Mc
% Mc1 = Mc1./Mc; Mc1 = Mc1/max(max(abs(Mc1))); % normalize values between -1 and 1
% Mab = Mab./Mc; Mab = Mab/max(max(abs(Mab))); % normalize values between -1 and 1
end