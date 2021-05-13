% ICSCompiler(...) is meant for fitting ICS autocorrelation functions
% (ACFs), and to return information about the fit(s).
%
% INPUT PARAMS
% J: the movie to be analyzed (1st 2 dims are spatial, 3rd is temporal)
% xi(eta)_lags: spatial lags to fit in x(y)-direction (note x is dim 2, and
%               y is dim 1).
%

function [ics_run] = ICSCompiler(J,xi_lags,eta_lags,varargin)
%% varargin

% boolean to include (0,0) spatial lag when fitting
includeZeroLag = 0;
% number of subsets to split ACF into (fit is run on each subset)
n_ics_subs = 1;
% mean type for subtraction in ICS  
% temporal mean subtraction is useful for spatially inhomogeneous data
mean_type = 'temporal';
% 
% background correction logical
bkgd_cxn = 0;
% background type (only relevant if bkgd_cxn=1)
bkgd_type = 'mean';
%
% path for saving best fit figures of ICS subsets
figpath = '';
% show best fit figures of ICS subsets; has no effect if figures are being
% saved
show_figs = 0;
for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'includeZeroLag','zeroLag'}))
        if varargin{ii+1} == 0 || varargin{ii+1} == 1
            includeZeroLag = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'nSubsets','subsets'}))
        if varargin{ii+1} >= 1
            n_ics_subs = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'meanType'}))
        if any(strcmpi(varargin{ii+1},{'spatial','temporal'}))
            mean_type = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'bkgdMovie','bkgdMov','noiseMovie',...
            'noiseMov'}))
        if isnumeric(varargin{ii+1}) 
            bkgd_cxn = 1;
            J_bkgd = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'bkgdType','noiseType'}))
        if any(strcmpi(varargin{ii+1},{'mean','bleach'}))
            bkgd_type = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'saveFigs'}))
        if ischar(varargin{ii+1})
            figpath = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'showFigs'}))
        show_figs = 1;
        if exist('varargin{ii+1}','var') && isnumeric(varargin{ii+1}) && varargin{ii+1} > 0
            n_figs = varargin{ii+1};
        end
    end
end

% default number of figures to show (if option selected)
if show_figs && ~exist('n_figs','var')
    n_figs = min(10,n_ics_subs); % max 10 figs if n_figs not specified
    fig_mod = round(n_ics_subs/n_figs); % modulo for figure showing
elseif exist('n_figs','var')    
    fig_mod = round(n_ics_subs/n_figs);
end
%

%% background noise

if bkgd_cxn
    bkgd_noise = getBkgd(J_bkgd,'method',bkgd_type);
    J_sub = J - bkgd_noise;
end

%% ics computation

% number of physical cores
num_cores = feature('numcores');

% total number of frames
T = size(J,3);
% time vector
t = 1:T;

% compute ICS ACF
if T == 1
    mean_type = 'spatial';
end

% number of frames in each subset (except maybe last one)
[sub_rng,n_ics_subs] = partitionArr(t,n_ics_subs);

% array to store acf subsets
corr_n = cell(1,n_ics_subs);
sub_corr_n = cell(1,n_ics_subs);
for ii = 1:n_ics_subs    
    % get time range for the iith ics subset
    sub_rng_ii = sub_rng{ii};
    
    % calculate ics acf over specified time range
    % note: not calculating ics acf in one pass is important if using
    % temporal mean subtraction to avoid subtracting temporal mean of
    % entire image series
    if bkgd_cxn % use background corrected normalization
        [corr_n{ii},xi_grid,eta_grid] = ICS(J(:,:,sub_rng_ii),...
            'meantype',mean_type,'bkgd',J_sub(:,:,sub_rng_ii));
    else % no background correction
        [corr_n{ii},xi_grid,eta_grid] = ICS(J(:,:,sub_rng_ii),...
            'meantype',mean_type);
    end
    % pick out a subset of the ics acf 
    [sub_corr_n{ii},xi_sub_grid,eta_sub_grid] = getICSCorrSub(corr_n{ii},xi_lags,...
        eta_lags,'includeZeroLag',includeZeroLag);
    if ii == 1
        % initialize ics acf subset average array 
        sub_corr_n_avg = zeros(size(sub_corr_n{ii},1),...
            size(sub_corr_n{ii},2),n_ics_subs);
    end
    % average ics acf subset
    sub_corr_n_avg(:,:,ii) = mean(sub_corr_n{ii},3);
end

%% fitting

% function handles of LSF error and fit functions
fit_fun = @(params) ICSFit(params,xi_sub_grid,eta_sub_grid);
err = @(params,y_data) ICSFit(params,xi_sub_grid,eta_sub_grid,'err',y_data);

% array to store best fit params
opt_params = zeros(n_ics_subs,3);
% opt_params = zeros(n_ics_subs,n_params);

% array to store objective function values
err_min = zeros(1,n_ics_subs);

% arrays to store initial guesses and bounds
guess_params = zeros(n_ics_subs,3);
lb = zeros(n_ics_subs,3);
ub = zeros(n_ics_subs,3);
% guess_params = zeros(n_ics_subs,n_params);
% lb = zeros(n_ics_subs,n_params);
% ub = zeros(n_ics_subs,n_params);

err_ii = cell(size(sub_corr_n_avg));
for ii = 1:n_ics_subs
    err_ii{ii} = @(params) err(params,sub_corr_n_avg(:,:,ii));
end

disp('fitting ICS subsets...')
tic
if num_cores > 1 && n_ics_subs >= num_cores
    parfor ii = 1:n_ics_subs
        % get initial guess
        corr_ii = corr_n{ii};
        sub_corr_ii = sub_corr_n{ii};
        
        [guess_params(ii,:),lb(ii,:),ub(ii,:)] = getICSGuess(corr_ii,sub_corr_ii);
        %
        
        opts = optimoptions(@fmincon,'Algorithm','interior-point');
        problem = createOptimProblem('fmincon','objective',...
            err_ii{ii},'x0',guess_params(ii,:),'lb',lb(ii,:),'ub',ub(ii,:),'options',opts);
        gs = GlobalSearch('Display','off'); % global search object
        [opt_params(ii,:),err_min(ii)] = run(gs,problem);
    end
    delete(gcp)
else
    for ii = 1:n_ics_subs
        % get initial guess
        corr_ii = corr_n{ii};
        sub_corr_ii = sub_corr_n{ii};
        
        [guess_params(ii,:),lb(ii,:),ub(ii,:)] = getICSGuess(corr_ii,sub_corr_ii);
        %
        
%         guess_params = [guess_params(1),guess_params(3),guess_params(2),guess_params(2),rand*pi];
%         lb=[lb(1),lb(3),lb(2),lb(2),eps];
%         ub=[ub(1),ub(3),ub(2),ub(2),pi];
%         
        opts = optimoptions(@fmincon,'Algorithm','interior-point');
        problem = createOptimProblem('fmincon','objective',...
            err_ii{ii},'x0',guess_params(ii,:),'lb',lb(ii,:),'ub',ub(ii,:),'options',opts);
        gs = GlobalSearch('Display','off'); % global search object
        [opt_params(ii,:),err_min(ii)] = run(gs,problem);
    end
end
toc

% save raw ICS autocorr data
ics_run.corr_n = corr_n; % ics_run.corr = corr;
ics_run.xi_grid = xi_grid; ics_run.eta_grid = eta_grid;
ics_run.sub_corr_n = sub_corr_n; % ics_run.sub_corr = sub_corr;
ics_run.xi_sub_grid = xi_sub_grid; ics_run.eta_sub_grid = eta_sub_grid;
ics_run.sub_corr_t = sub_corr_n_avg;

% save guess params and bounds
ics_run.guess_params = guess_params;
ics_run.lb = lb;
ics_run.ub = ub;

% save fit results
ics_run.opt_params = opt_params;
ics_run.err_min = err_min;

% save miscellaneous
ics_run.n_ics_subs = n_ics_subs;

%% plot and save figures

if show_figs && ~exist('n_figs','var')
    % n_figs unspecified; show max 10 figs
    n_figs = min(10,n_ics_subs);
end
% fig_inds is a partitioned enumeration of all figures; first figure
% index in each partition is displayed
[fig_inds,n_figs] = partitionArr(n_ics_subs,n_figs); 
%

% check if correlation was queried as vector or array
xi_single = isscalar(xi_lags); eta_single = isscalar(eta_lags);
if any([xi_single,eta_single])
	is_vec = 1;
else
    is_vec = 0;
end
%

if ~isempty(figpath) || show_figs
    if ~isempty(figpath) && ~exist(figpath,'dir')
        mkdir(figpath)
    end
    
    [y_mid,x_mid] = getCtrPxl(corr_n{1});
    [y_mid_sub,x_mid_sub] = getCtrPxl(sub_corr_n{1});
    
    xi_vec = reshape(xi_sub_grid,[1,numel(xi_sub_grid)]);
    eta_vec = reshape(eta_sub_grid,[1,numel(eta_sub_grid)]);
    for ii = 1:n_figs
        show_ind = fig_inds{ii}(1);
        zero_lag_avg_n = mean(corr_n{show_ind}(y_mid,x_mid,:));
        
        figure()
        hold on       
        if ~is_vec
            corr_vec = reshape(sub_corr_n_avg(:,:,show_ind),[1,numel(sub_corr_n_avg(:,:,show_ind))]);
            scatter3(xi_vec,eta_vec,corr_vec,'ok')
            scatter3(0,0,zero_lag_avg_n,'ok')
            surf(xi_sub_grid,eta_sub_grid,fit_fun(opt_params(show_ind,:)))
            %         mesh(xi_sub_grid,eta_sub_grid,fit_fun(guess_params(ii,:)))
        else
            plot(xi_lags+eta_lags,sub_corr_n_avg(y_mid_sub+eta_lags,x_mid_sub+xi_lags,show_ind),...
                '.','markersize',16)
            plot(0,zero_lag_avg_n,'.','markersize',16)
            plot(xi_lags+eta_lags,ICSFit(opt_params(show_ind,:),xi_lags,eta_lags),'-','linewidth',2)
            plot(xi_lags+eta_lags,ICSFit(guess_params(show_ind,:),xi_lags,eta_lags),'--','linewidth',2)
        end

        xlabel('$\xi$ (pixels)','fontsize',12,'interpreter','latex')
        ylabel('$\phi(\xi)$','fontsize',12,'interpreter','latex')
        legend({'data','zero lag','fit','guess'},'fontsize',12,'interpreter','latex')
        
        tightfig(gcf)
        
        if ~isempty(figpath)
            filename = [figpath,filesep,'ics_sub_fit_',num2str(show_ind),'.fig'];
            saveas(gcf,filename);
            filename = [figpath,filesep,'ics_sub_fit_',num2str(show_ind),'.pdf'];
            saveas(gcf,filename);
            
            close(gcf)
        end
    end
end