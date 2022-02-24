function SDL_WiringCost_Th_Age10(SDL)

%==========================================================================
%   (1) calculate minimum wiring cost, and
%   (2) calculate the binary mediation matrix after thresholding
% Input
% ------ SDL, all important parameters
% Output
% ------ SDL, including updated
%        SDL.WC, minimum wiring cost of each matrix (row-group, column-C, C1 & ab)
%        SDL.Thr, threshold corresponding to wc

%% Load data
fdir = fullfile(SDL.out,[SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4}]);
fn = fullfile(fdir,['Matrix_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Matrix in %s\n',fn);

W = []; % to contain the wiring cost
%% minimum wiring cost & threshold
if strcmp(SDL.ana_type{1},'med') || strcmp(SDL.ana_type{1},'medcov')% mediation matrix
    k = 0.5 * size(MMc,1);
    atype = 'dir'; % connection with direction
    % abslute values
    Thr = []; WC = []; ftype = 'abs';
    D = MMc(1:k,:,1);      [Thr(1,1), WC(1,1)] = SDL_ThresMinEdges(D,ftype,atype); % for M1c
    D = MMc(k+1:end,:,1);  [Thr(2,1), WC(2,1)] = SDL_ThresMinEdges(D,ftype,atype); % for M2c
    D = MMc1(1:k,:,1);     [Thr(1,2), WC(1,2)] = SDL_ThresMinEdges(D,ftype,atype); % for M1c1
    D = MMc1(k+1:end,:,1); [Thr(2,2), WC(2,2)] = SDL_ThresMinEdges(D,ftype,atype); % for M2c1
    D = MMab(1:k,:,1);     [Thr(1,3), WC(1,3)] = SDL_ThresMinEdges(D,ftype,atype); % for M1ab
    D = MMab(k+1:end,:,1); [Thr(2,3), WC(2,3)] = SDL_ThresMinEdges(D,ftype,atype); % for M2ab
    W.WC_abs = WC; W.Thr_abs = Thr;
    
    % positive values
    Thr = []; WC = []; ftype = 'pos';
    D = MMc(1:k,:,1);      [Thr(1,1), WC(1,1)] = SDL_ThresMinEdges(D,ftype,atype); % for M1c
    D = MMc(k+1:end,:,1);  [Thr(2,1), WC(2,1)] = SDL_ThresMinEdges(D,ftype,atype); % for M2c
    D = MMc1(1:k,:,1);     [Thr(1,2), WC(1,2)] = SDL_ThresMinEdges(D,ftype,atype); % for M1c1
    D = MMc1(k+1:end,:,1); [Thr(2,2), WC(2,2)] = SDL_ThresMinEdges(D,ftype,atype); % for M2c1
    D = MMab(1:k,:,1);     [Thr(1,3), WC(1,3)] = SDL_ThresMinEdges(D,ftype,atype); % for M1ab
    D = MMab(k+1:end,:,1); [Thr(2,3), WC(2,3)] = SDL_ThresMinEdges(D,ftype,atype); % for M2ab
    W.WC_pos = WC; W.Thr_pos = Thr;
    
    % negative values
    Thr = []; WC = []; ftype = 'neg';
    D = MMc(1:k,:,1);      [Thr(1,1), WC(1,1)] = SDL_ThresMinEdges(D,ftype,atype); % for M1c
    D = MMc(k+1:end,:,1);  [Thr(2,1), WC(2,1)] = SDL_ThresMinEdges(D,ftype,atype); % for M2c
    D = MMc1(1:k,:,1);     [Thr(1,2), WC(1,2)] = SDL_ThresMinEdges(D,ftype,atype); % for M1c1
    D = MMc1(k+1:end,:,1); [Thr(2,2), WC(2,2)] = SDL_ThresMinEdges(D,ftype,atype); % for M2c1
    D = MMab(1:k,:,1);     [Thr(1,3), WC(1,3)] = SDL_ThresMinEdges(D,ftype,atype); % for M1ab
    D = MMab(k+1:end,:,1); [Thr(2,3), WC(2,3)] = SDL_ThresMinEdges(D,ftype,atype); % for M2ab
    W.WC_neg = WC; W.Thr_neg = -1*Thr;
    
else % correlation or regression matrix
    if strcmp(SDL.ana_type{1},'corr') ||strcmp(SDL.ana_type{1},'partialcorr')  % correlation matrix
        MMc = MRHO;
    elseif strcmp(SDL.ana_type{1},'regr') || strcmp(SDL.ana_type{1},'regrcov') % regression matrix
        MMc = MB;
    else
    end
    atype = 'und'; % connection without direction
    k = size(MMc,1)/16; % each group has 148 regions
    % abslute values
    Thr = []; WC = []; ftype = 'abs';
    D = MMc(    1:k,  :,1);  [Thr(1,1), WC(1,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g11
    D = MMc(  k+1:2*k,:,1);  [Thr(2,1), WC(2,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g12
    D = MMc(2*k+1:3*k,:,1);  [Thr(3,1), WC(3,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g13
    D = MMc(3*k+1:4*k,:,1);  [Thr(4,1), WC(4,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g14
    D = MMc(4*k+1:5*k,:,1);  [Thr(5,1), WC(5,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g15
    D = MMc(5*k+1:6*k,:,1);  [Thr(6,1), WC(6,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g16
    D = MMc(6*k+1:7*k,:,1);  [Thr(7,1), WC(7,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g17
    D = MMc(7*k+1:8*k,:,1);  [Thr(8,1), WC(8,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g18
    
    D = MMc(8*k+1:9*k,:,1);  [Thr(9,1), WC(9,1)]  = SDL_ThresMinEdges(D,ftype,atype); % for g21
    D = MMc(9*k+1:10*k,:,1); [Thr(10,1),WC(10,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g22
    D = MMc(10*k+1:11*k,:,1);[Thr(11,1),WC(11,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g23
    D = MMc(11*k+1:12*k,:,1);[Thr(12,1),WC(12,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g24
    D = MMc(12*k+1:13*k,:,1);[Thr(13,1),WC(13,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g25
    D = MMc(13*k+1:14*k,:,1);[Thr(14,1),WC(14,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g26
    D = MMc(14*k+1:15*k,:,1);[Thr(15,1),WC(15,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g27
    D = MMc(15*k+1:16*k,:,1);[Thr(16,1),WC(16,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g28
    W.WC_abs = WC; W.Thr_abs = Thr;
    
    % positive values
    Thr = []; WC = []; ftype = 'pos';
    D = MMc(    1:k,  :,1);  [Thr(1,1), WC(1,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g11
    D = MMc(  k+1:2*k,:,1);  [Thr(2,1), WC(2,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g12
    D = MMc(2*k+1:3*k,:,1);  [Thr(3,1), WC(3,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g13
    D = MMc(3*k+1:4*k,:,1);  [Thr(4,1), WC(4,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g14
    D = MMc(4*k+1:5*k,:,1);  [Thr(5,1), WC(5,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g15
    D = MMc(5*k+1:6*k,:,1);  [Thr(6,1), WC(6,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g15
    D = MMc(6*k+1:7*k,:,1);  [Thr(7,1), WC(7,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g15
    D = MMc(7*k+1:8*k,:,1);  [Thr(8,1), WC(8,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g15
    
    D = MMc(8*k+1:9*k,:,1);    [Thr(9,1),  WC(9,1)]  = SDL_ThresMinEdges(D,ftype,atype); % for g21
    D = MMc(9*k+1:10*k,:,1);   [Thr(10,1), WC(10,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g22
    D = MMc(10*k+1:11*k,:,1);  [Thr(11,1), WC(11,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g23
    D = MMc(11*k+1:12*k,:,1);  [Thr(12,1), WC(12,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g24
    D = MMc(12*k+1:13*k,:,1);  [Thr(13,1), WC(13,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g25
    D = MMc(13*k+1:14*k,:,1);  [Thr(14,1), WC(14,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g26
    D = MMc(14*k+1:15*k,:,1);  [Thr(15,1), WC(15,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g27
    D = MMc(15*k+1:16*k,:,1);  [Thr(16,1), WC(16,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g28
    W.WC_pos = WC; W.Thr_pos = Thr;
    
    % negative values
    Thr = []; WC = []; ftype = 'neg';
    D = MMc(    1:k,  :,1);  [Thr(1,1), WC(1,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g11
    D = MMc(  k+1:2*k,:,1);  [Thr(2,1), WC(2,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g12
    D = MMc(2*k+1:3*k,:,1);  [Thr(3,1), WC(3,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g13
    D = MMc(3*k+1:4*k,:,1);  [Thr(4,1), WC(4,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g14
    D = MMc(4*k+1:5*k,:,1);  [Thr(5,1), WC(5,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g15
    D = MMc(5*k+1:6*k,:,1);  [Thr(6,1), WC(6,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g16
    D = MMc(6*k+1:7*k,:,1);  [Thr(7,1), WC(7,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g17
    D = MMc(7*k+1:8*k,:,1);  [Thr(8,1), WC(8,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g18
    
    D = MMc(8*k+1:9*k,:,1);    [Thr(9,1), WC(9,1)]   = SDL_ThresMinEdges(D,ftype,atype); % for g21
    D = MMc(9*k+1:10*k,:,1);   [Thr(10,1), WC(10,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g22
    D = MMc(10*k+1:11*k,:,1);  [Thr(11,1), WC(11,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g23
    D = MMc(11*k+1:12*k,:,1);  [Thr(12,1), WC(12,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g24
    D = MMc(12*k+1:13*k,:,1);  [Thr(13,1), WC(13,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g25
    D = MMc(13*k+1:14*k,:,1);  [Thr(14,1), WC(14,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g26
    D = MMc(14*k+1:15*k,:,1);  [Thr(15,1), WC(15,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g27
    D = MMc(15*k+1:16*k,:,1);  [Thr(16,1), WC(16,1)] = SDL_ThresMinEdges(D,ftype,atype); % for g28
    W.WC_neg = WC; W.Thr_neg = -1*Thr;
    
end
fn = fullfile(fdir,['WC_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
save(fn,'W'); fprintf('Saved: Minimum Wiring Cost & Corresponding Thresholds in %s\n',fn);

% end
end

%calculate the minimum wiring cost of corresponding threshold
function [Thr,WC] = SDL_ThresMinEdges(D,ftype,atype)
% Input
% ------ D, mediation matrix
% ------ ftype, 'pos', 'neg', 'abs'
% ------ atype, 'und','dir'
% Output
% ------ WC,  minimum wiring cost
% ------ Thr, the threshold corresponding to the minimum wiring cost
% Setting values below th into 0
if strcmp(ftype,'abs')
    D = abs(D); % absolute values
elseif strcmp(ftype,'pos')
    D(D<0) = 0; % keep positive values
elseif strcmp(ftype,'neg')
    D(D>0) = 0; % keep negative values
    D = -1*D; % for the ease of calculating
else
end
D1 = D;
D(1:length(D)+1:end)=0;% remove the values in diagonal
D1(1:length(D)+1:end)=NaN;% remove the values in diagonal
th = 0.5 * (min(min(D1)) + max(max(D1))); % threshold changes beginning from the middle
Lstep = 0.0001 * (max(max(D1)) - min(min(D1))); % delta per step
N = 10000; % maximum steps
j = 0;
flag = 0; % 0 if the threshold has not been detected, 1 if it is calculated
while (j < N) && (flag == 0)
    j = j + 1;
    Du = D; Dd = D; % two matrixes corresponding to higher and lower thresholds
    Du(Du<th+Lstep)=0; Dd(Dd<th-Lstep)=0; % values below threshold are set into 0
    Du(Du>=th+Lstep)=1;Dd(Dd>=th-Lstep)=1;% values above or equal to threshold are set into 1
    if strcmp(atype,'und') % without direction, e.g. corr or regr
        deg_Du = degrees_und(Du);
        deg_Dd = degrees_und(Dd);
    else % with direction, e.g. mediation
        [id,od,deg_Du] = degrees_dir(Du);
        [id,od,deg_Dd] = degrees_dir(Dd);
    end
    %     if (sum(any(Du)) < length(Du)) && (sum(any(Dd)) < length(Dd))
    if   ~isempty(find(~deg_Du)) && ~isempty(find(~deg_Dd))
        % both up and down matrix have isolated nodes,
        % i.e. too strict, need to lower thresholds
        th = th - Lstep; % reduce threshold if too sparse
        flag = 0; % need further search
    elseif isempty(find(~deg_Du)) && isempty(find(~deg_Dd))
        % both up and down matrix do not have isolated nodes,
        % i.e. too loose, and need to increase thresholds
        th = th + Lstep; % increase threshold if redundant connections
        flag = 0; % need further search
    else % the threshold is detected
        flag = 1;
    end
end

%     %% calculating the threshold, i.e. th, through Prof. He and Dr. Wang's way
%     D1 = D; % e.g., AllCorrMat is the partial correlation matrix
%     D1(D < Thr) = 0; % Thr is a threshold, e.g. 0.20
%     [ci, sizes] = components(sparse(D1));
%     % should plot the relationship between Thr and maximum numbers of components

% report threshold
if flag == 0 % fail to find the threshold
    Thr = Inf;
else
    Thr = th;
end

Data  = double(D>Thr);

if strcmp(atype,'und') % without direction, e.g. corr or regr
    [kden,N,K] = density_und(Data); 
else % with direction, e.g. mediation
    [kden,N,K] = density_dir(Data); 
end
WC = kden;

end