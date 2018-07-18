%%% Initialize program %%%

disp(' ')
disp('%%%% MCMC Santa Cruz Mountains Single K %%%%')
disp(' ')
t_tot=0;
c=clock;
PID=feature('getpid');
disp(['Date/time in Stanford= ',num2str(c(1)),' ',num2str(c(2)),' ',...
    num2str(c(3)),' ',num2str(c(4)),' ',num2str(c(5)),' ',num2str(c(6))])
disp(['PID= ',num2str(PID)])

rng(12309);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load up data: %%%%

tic
% Change working directory accordingly
addpath(genpath('/data/cees/aron/SCM/MCMC'))
% addpath(genpath('/data/cees/aron/bin/matlab'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Comment for server runs %%%%%
% if script with input variables %
%%%% is available
%%% Uncomment for local testing %%%

% v=0.25;         %Poisson's ratio
% charlen=0.25;   %Mesh characteristic length
% theta = 0.4;    %Concavity (slope of log stream power law profiles)
% 
% sig_elevation = 60;
% sig_bay= 5e-5;
% 
% Niter = 1e4;
% 
% % Parse out model initiual guess, step size, bounds and prior
% X.x0     = [0; 0; -4];
% X.xstep  = 0.1*[1e-5; 1e-5; 1e-3];
% X.xbnds  = [-1e-1 1e-1; -1e-1 1e-1; -8 -2];
% X.xprior = [];
% X.C = Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp_K_unit=1-2*theta;

disp(['Poisson''s ratio tectonic model= ',num2str(v)])
disp(['Mesh characteristic length tectonic model= ',num2str(charlen)])
disp(['Concavity stream profiles (theta)= ',num2str(theta)])
disp('MCMC model parameters:')
disp(['- initial guess= ',num2str(X.x0(1)),'; ',num2str(X.x0(2)),...
    '; ',num2str(X.x0(3))])
disp(['   shear; push; log10(K) [m/yr; m/yr; log10(m^',num2str(exp_K_unit),'/yr)]'])
disp(['- step size= ',num2str(X.xstep(1)),'; ',num2str(X.xstep(2)),...
    '; ',num2str(X.xstep(3))])
disp(['   shear; push; log10(K) [m/yr; m/yr; log10(m^',num2str(exp_K_unit),'/yr)]'])
disp(['- bounds= ',num2str(X.xbnds(1,1)),' ',num2str(X.xbnds(1,2)),'; ',...
    num2str(X.xbnds(2,1)),' ',num2str(X.xbnds(2,2)),'; ',...
    num2str(X.xbnds(3,1)),' ',num2str(X.xbnds(3,2))])
disp(['   shear; push; log10(K) [m/yr; m/yr; log10(m^',num2str(exp_K_unit),'/yr)]'])
disp(['- prior= ',num2str(X.xprior)])
disp(['- C= ',num2str(X.C)])
disp(['Sigma elevation= ',num2str(sig_elevation),' m'])
disp(['Sigma Bay area 0 uplift constrain= ',num2str(sig_bay),' mm/yr'])
disp(['Number of MCMC realizations= ',num2str(Niter)])

[geo_map, d, Ginv_elev, channel_indexes, channel_elevations, n_K, compare_elements, litho_chan] = prepare_inputs_for_objective_function(theta);
channel_elevations_to_compare = channel_elevations(compare_elements);

disp(['Number of Geo units= ',num2str(n_K)])

G_tect=load('G_tect_shear_push.txt');
G_bay=load('G_bay_shear_push.txt');

disp(['Number of channel observation points= ',num2str(length(G_tect(:,1)))])

% Take account of the 2 BC due to 1 imposed on each side of the box
G_tect=G_tect./2; G_bay=G_bay./2;

%%%%%%% Parse out data and errors

% D.d = channel_elevations_to_compare;
% Adding bay area 0 uplift-rate
D.d = [channel_elevations_to_compare;zeros(length(G_bay(:,1)),1)];

Sig_vec=[(sig_elevation^2)*ones(length(channel_elevations_to_compare),1);...
    (sig_bay^2)*ones(length(G_bay(:,1)),1)];

D.Sig = sparse(diag(Sig_vec));

% Parse out model initiual guess, step size, bounds and prior

% X.x0     = [0; 0; -4];
% X.xstep  = 0.1*[1e-5; 1e-5; 1e-3];
% X.xbnds  = [-1e-1 1e-1; -1e-1 1e-1; -8 -2];
% X.xprior = [];
% X.C = Inf;

% Niter = 2e6;

t=toc;
t_tot=t_tot+t;
disp(['Loading files, parsing out data, errors, and model takes ',num2str(t),' secs'])


%%%%% Run MCMC routine %%%%%

tic

[x_keep, logL_keep, logQ_keep, accrate] = mcmc_new('elev_fun_ibem_n1_platemotion_Bay_1K',D,X,Niter,1,G_tect,G_bay,d,Ginv_elev,channel_indexes,compare_elements);

t=toc;
t_tot=t_tot+t;
disp(['Running MCMC routine takes ',num2str(t),' secs'])


%%%%% Save MCMC %%%%%

tic
filename = ['mcmc_scm_singleK_',num2str(c(1)),'_',num2str(c(2)),'_',...
    num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)),'_',num2str(c(6))];
% filename = strcat(filename(1:end-2),['_PST_PID',num2str(PID),'_v',num2str(v),...
%     '_cl',num2str(charlen),'_theta',num2str(theta),'_SigElev',num2str(sig_elevation),...
%     '_SigBay',num2str(sig_bay)]);
filename = strcat(filename(1:end-2),['_PST_PID',num2str(PID),'.mat']);

save(filename, 'x_keep', 'logL_keep', 'logQ_keep', 'accrate', '-v7.3')
t=toc;
t_tot=t_tot+t;
disp(['Saving files takes ',num2str(t),' secs'])

disp(['Program takes a total of ',num2str(t_tot),' secs'])
disp('good bye!')

quit