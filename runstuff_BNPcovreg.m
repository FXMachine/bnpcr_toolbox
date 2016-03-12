%% SCRIPT FOR RUNNING BASIC SAMPLER ON SIMULATED DATA EXAMPLE %%

% Load simulated data:
load simData.mat

% y is a p x N matrix specifying the N observations of the p-dimensional data
[p N] = size(y);

% Indicate which observations are present.  A 0 indicates a missing value.
% Note that the following script will impute any missing values.  For
% analytic marginalization of missing values, run BNP_covreg_varinds.
% See runstuff_varinds_flu.m for an example in this case.
inds_y = ones(size(y));
inds_y = inds_y > 0;

% Rescale predictor space for forming correlation matrix:
x = [1:N]./N;

c = 10;  % Gaussian process length scale parameter.  Smaller c will lead to slower changes in correlations.
d = 1;
r = 1e-5;
K = zeros(N);
for ii=1:N
    for jj=1:N
        dist_ii_jj = abs(x(ii)-x(jj));
        K(ii,jj) = d*exp(-c*(dist_ii_jj^2));
    end
end
K = K + diag(r*ones(1,N));
invK = inv(K);
logdetK = 2*sum(log(diag(chol(K))));

prior_params.K.c_prior = 1;
prior_params.K.invK = invK;
prior_params.K.K = K;
prior_params.K.logdetK = logdetK;
prior_params.sig.a_sig = 1;
prior_params.sig.b_sig = 0.1;
prior_params.hypers.a_phi = 1.5;
prior_params.hypers.b_phi = 1.5;
prior_params.hypers.a1 = 2; 
prior_params.hypers.a2 = 2; 

settings.L = 5;    % truncation level for dictionary of latent GPs
settings.k = 4;    % latent factor dimension
settings.Niter = 10000;  % number of Gibbs iterations to run
settings.saveEvery = 100;  % how often to write stats to disk
settings.storeEvery = 10;  % how often to store samples in struct
settings.saveMin = 1;
settings.saveDir = 'test';  % name of directory that will be created to save stats in
settings.trial = 1;  % trial number specified here provides an index for labeling saved files
settings.sample_K_flag = 3; % choose to fix bandwidth parameter instead of sampling it
settings.latent_mean = 1; % indicate whether one wants to model a non-zero latent mean \mu(x)
settings.inds_y = inds_y;

% Indicate whether we are restarting from a previous chain:
restart = 0; 
if restart
    settings.lastIter = 10000; % last saved iteration we want to load
end

BNP_covreg_varinds(y,prior_params,settings,restart,true_params);

%%

% If one wants to sample the bandwidth parameter, specify some vector of
% possible values in c_vec and use the following in place of the
% corresponding code above.

c = 10;  % Gaussian process length scale parameter.  Smaller c will lead to slower changes in correlations.
d = 1;
r = 1e-5;
c_vec = c;
Nc = length(c_vec);
invK = zeros(N,N,Nc);
K = zeros(N,N,Nc);
logdetK = zeros(1,Nc);
for kk=1:Nc
    K_tmp = zeros(N,N);
    for ii=1:N
        for jj=1:N
            dist_ii_jj = abs(x(ii)-x(jj));
            K_tmp(ii,jj) = d*exp(-c_vec(kk)*(dist_ii_jj^2));
        end
    end
    K_tmp = K_tmp + diag(r*ones(1,N));
    K(:,:,kk) = K_tmp;
    invK(:,:,kk) = inv(K_tmp);
    logdetK(kk) = 2*sum(log(diag(chol(K_tmp))));
end

prior_params.K.c_prior = ones(1,Nc)/Nc;
prior_params.K.invK = invK;
prior_params.K.K = K;
prior_params.K.logdetK = logdetK;

% There are two options for how to sample K...

% (1) Sample K based on the p x N dimensional Gaussian posterior
% p(K | y, theta, eta, Sigma_0) having marginalized the latent GP
% functions zeta:
settings.sample_K_flag = 1; 

% (2) If p x N is too large, the above is infeasible.  Instead, sample the
% GP cov matrix K given the latent GP functions zeta p(K | zeta) (yielding 
% K cond ind of everything else, though leading to much slower mixing rates)
settings.sample_K_flag = 2; 
