function SDL_matrix_Age10(SDL)

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

idx11 = strcmp(T.Group,'PTSD') & (T.Age<10);              [idx11,~] = find(idx11==1); % group 11, PTSD & Age < 10
idx12 = strcmp(T.Group,'PTSD') & (T.Age>=10 & T.Age<=21);  [idx12,~] = find(idx12==1);% group 12, PTSD & 10<=Age<=21
idx13 = strcmp(T.Group,'PTSD') & (T.Age>=22 & T.Age<=39);  [idx13,~] = find(idx13==1);% group 13, PTSD & 22<=Age<=39
idx14 = strcmp(T.Group,'PTSD') & (T.Age>=40 & T.Age<=59);  [idx14,~] = find(idx14==1);% group 14, PTSD & 40<=Age<=59
%idx15 = strcmp(T.Group,'PTSD') & (T.Age>=60 & T.Age<40);  [idx15,~] = find(idx15==1);% group 15, PTSD & 30<=Age<40
%idx16 = strcmp(T.Group,'PTSD') & (T.Age>=40 & T.Age<50);  [idx16,~] = find(idx16==1);% group 16, PTSD & 40<=Age<50
%idx17 = strcmp(T.Group,'PTSD') & (T.Age>=50 & T.Age<60);  [idx17,~] = find(idx17==1);% group 17, PTSD & 50<=Age<60
idx18 = strcmp(T.Group,'PTSD') & (T.Age>=60);             [idx18,~] = find(idx18==1);% group 18, PTSD & 60<=Age

idx21 = strcmp(T.Group,'CONT') & (T.Age<10);              [idx21,~] = find(idx21==1);% group 21, CONT & Age < 10
idx22 = strcmp(T.Group,'CONT') & (T.Age>=10 & T.Age<=21);  [idx22,~] = find(idx22==1);% group 22, CONT & 10<=Age<=21
idx23 = strcmp(T.Group,'CONT') & (T.Age>=22 & T.Age<=39);  [idx23,~] = find(idx23==1);% group 23, CONT & 22<=Age<=39
idx24 = strcmp(T.Group,'CONT') & (T.Age>=40 & T.Age<=59);  [idx24,~] = find(idx24==1);% group 24, CONT & 40<=Age<=59
%idx25 = strcmp(T.Group,'CONT') & (T.Age>=30 & T.Age<40);  [idx25,~] = find(idx25==1);% group 25, CONT & 30<=Age<40
%idx26 = strcmp(T.Group,'CONT') & (T.Age>=40 & T.Age<50);  [idx26,~] = find(idx26==1);% group 26, CONT & 40<=Age<50
%idx27 = strcmp(T.Group,'CONT') & (T.Age>=50 & T.Age<60);  [idx27,~] = find(idx27==1);% group 27, CONT & 50<=Age<60
idx28 = strcmp(T.Group,'CONT') & (T.Age>=60);             [idx28,~] = find(idx28==1);% group 28, CONT & 60<=Age

[   length(idx11),length(idx21);...
    length(idx12),length(idx22);...
    length(idx13),length(idx23);...
    length(idx14),length(idx24);...
    %length(idx15),length(idx25);...
    %length(idx16),length(idx26);...
    %length(idx17),length(idx27);...
    length(idx18),length(idx28)]


%% (2) permutation N times by shuffling group labels

N = SDL.N; % times of permutation: N>1, calculate the raw mediation matrix,
%           as well as N-1 mediation matrix with shuffling labels between
%           groups

% if strcmp(SDL.ana_type{1},'corr') || strcmp(SDL.ana_type{1},'partialcorr') % correlation/partial-correlation matrix
    if strcmp(SDL.ana_type{1},'corr')
        SDL_corr = @corr;        % linear partial correlation coefficients between pairs of variables without controlling other variables
    else
        SDL_corr = @partialcorr; % linear partial correlation coefficients between pairs of variables with controlling other variables
    end
    for i = 1:N+1 % per permutation
        tic;
        if i == 1 % no permutation
            RHO11 = SDL_corr(R(idx11,:)); % linear partial correlation coefficients between pairs of variables in group 1
            RHO12 = SDL_corr(R(idx12,:)); % linear partial correlation coefficients between pairs of variables in group 2
            RHO13 = SDL_corr(R(idx13,:));
            RHO14 = SDL_corr(R(idx14,:));
            %RHO15 = SDL_corr(R(idx15,:));
            %RHO16 = SDL_corr(R(idx16,:));
            %RHO17 = SDL_corr(R(idx17,:));
            RHO18 = SDL_corr(R(idx18,:));
            
            RHO21 = SDL_corr(R(idx21,:)); % linear partial correlation coefficients between pairs of variables in group 1
            RHO22 = SDL_corr(R(idx22,:)); % linear partial correlation coefficients between pairs of variables in group 2
            RHO23 = SDL_corr(R(idx23,:));
            RHO24 = SDL_corr(R(idx24,:));
            %RHO25 = SDL_corr(R(idx25,:));
            %RHO26 = SDL_corr(R(idx26,:));
            %RHO27 = SDL_corr(R(idx27,:));
            RHO28 = SDL_corr(R(idx28,:));
        else % permutation by shuffling group labels
            [RHO11,RHO21] = SDL_rand(R,SDL_corr,idx11,idx21);
            [RHO12,RHO22] = SDL_rand(R,SDL_corr,idx12,idx22);
            [RHO13,RHO23] = SDL_rand(R,SDL_corr,idx13,idx23);
            [RHO14,RHO24] = SDL_rand(R,SDL_corr,idx14,idx24);
            %[RHO15,RHO25] = SDL_rand(R,SDL_corr,idx15,idx25);
            %[RHO16,RHO26] = SDL_rand(R,SDL_corr,idx16,idx26);
            %[RHO17,RHO27] = SDL_rand(R,SDL_corr,idx17,idx27);
            [RHO18,RHO28] = SDL_rand(R,SDL_corr,idx18,idx28);
        end
        MRHO(:,:,i) = [RHO11;RHO12;RHO13;RHO14;RHO18;...
                       RHO21;RHO22;RHO23;RHO24;RHO28];
        fprintf('Completed: permutation %d/%d for correlation matrix\t',i,N); toc;
    end
    fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
    fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
%     SDL.CorrMatrix = fn; % store the name of Results_MediationMatrix_xxx.mat file
    save(fn,'-v7.3'); fprintf('Saved: Matrix in %s\n',fdir);
    
% elseif strcmp(SDL.ana_type{1},'regr') % regression matrix
%     for i = 1:N % per permutation
%         tic;
%         if i == 1 % no permutation
%             B1 = SDL_regression(R(g1,:)); % regression coefficients between pairs of variables in group 1
%             B2 = SDL_regression(R(g2,:)); % regression coefficients between pairs of variables in group 1
%         else % permutation by shuffling group labels
%             rmlist = randperm(length(g1)+length(g2));
%             B1 = SDL_regression(R(rmlist(g1),:)); % regression coefficients between pairs of variables in group 1
%             B2 = SDL_regression(R(rmlist(g2),:)); % regression coefficients between pairs of variables in group 1
%         end
% %         MB(:,:,i) = [B1;B2]; 
%         MB(:,:,i) = [(B1+B1')/2;(B2+B2')/2]; % to make the matrix symmetric
%         fprintf('Completed: permutation %d/%d\t',i,N); toc;
%     end
%     fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
%     fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
%     SDL.RegrMatrix = fn; % store the name of Results_MediationMatrix_xxx.mat file
%     save(fn,'-v7.3'); fprintf('Saved: Regression Matrix in %s\n',fn);
%     
% elseif strcmp(SDL.ana_type{1},'med') % mediation matrix
%     for i = 1:N % per permutation
%         tic;
%         if i == 1 % no permutation
%             [M1c,M1c1,M1ab] = SDL_mdiation(R(g1,:),RMV(g1,:)); % mediation matrix for group 1
%             [M2c,M2c1,M2ab] = SDL_mdiation(R(g2,:),RMV(g2,:)); % mediation matrix for group 2
%         else % permutation by shuffling group labels
%             rmlist = randperm(length(g1)+length(g2));
%             [M1c,M1c1,M1ab] = SDL_mdiation(R(rmlist(g1),:),RMV(rmlist(g1),:)); % mediation matrix for group 1
%             [M2c,M2c1,M2ab] = SDL_mdiation(R(rmlist(g2),:),RMV(rmlist(g2),:)); % mediation matrix for group 2
%         end
%         MMc(:,:,i) = [M1c;M2c]; MMc1(:,:,i) = [M1c1;M2c1]; MMab(:,:,i) = [M1ab;M2ab];
%         fprintf('Completed: permutation %d/%d\t',i,N); toc;
%     end
%     fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
%     fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
%     SDL.MedMatrix = fn; % store the name of Results_MediationMatrix_xxx.mat file
%     save(fn,'-v7.3'); fprintf('Saved: Mediation Matrix in %s\n',fn);
%     
% else
% end
disp('Matrix OK!');
% end
%end


% permutation by shuffling the group labels of two groups
% and then make corelation or partia correlation matrix
function [R1,R2] = SDL_rand(R,fun,idx1,idx2)
% Input
% ------ R, Matrix containing the raw data
% ------ fun, @corr or @partialcorr
% ------ idx1 & idx2, the index of subjects with each of the two groups
% Output
% ------ R1 & R2, correlation or partial correlation matrix of the permuted two groups
list0 = [idx1;idx2];
rl = randperm(length(list0));
list1 = list0(rl);% shuffle group labels
idx1a = list1(1:length(idx1));
idx2a = list1(1+length(idx1):end);
R1 = fun(R(idx1a,:));
R2 = fun(R(idx2a,:));
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
%             Cov = D; Cov(:,[i,j]) = []; % Covariates are too time conlengthing (~200 sec per permutation)
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
%             Cov = D; Cov(:,[i,j]) = []; % Covariates are too time conlengthing (~200 sec per permutation)
            B = regress(Y,[ones(size(X)),X]); 
            Mc(i,j)  = B(2);
            B = regress(Y,[ones(size(X)),X,M]); 
            Mc1(i,j) = B(2); 
            Mab(i,j) = Mc(i,j) - Mc1(i,j);
        end
    end
end
end
% % Below, Mc1 and Mab are not raw values but proportion of Mc
% Mc1 = Mc1./Mc; Mc1 = Mc1/max(max(abs(Mc1))); % normalize values between -1 and 1
% Mab = Mab./Mc; Mab = Mab/max(max(abs(Mab))); % normalize values between -1 and 1

end