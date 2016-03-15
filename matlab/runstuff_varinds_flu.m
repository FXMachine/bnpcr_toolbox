%%
    
% Script to run the sampler on the flu data.  Also provides an example of
% the necessary inputs for running the code that analytically marginalizes
% missing data.

load flu_US;

month_names = {'January','February','March','April','May','June','July','August','September','October','November','December'};

times = datenum(dates);
[years months days] = datevec(dates);

flu = data';

start_dates = zeros(1,size(flu,1));
for ii=1:size(flu,1)
    start_date_ii = 1;
    for tt=1:size(flu,2)
        if isnan(flu(ii,tt))
            start_date_ii = start_date_ii + 1;
        end
    end
    start_dates(ii) = start_date_ii;
end

[q T] = size(flu);

vars = zeros(q,1);
for i=1:q
    vars(i) = var(flu(i,start_dates(i):end));
end

flu = flu./sqrt(max(vars));

y=zeros(q,T);
for i=1:q
    y(i,start_dates(i):end) = flu(i,start_dates(i):end);
end

y = 1.75*y;

tmp = cumsum(sum(y,1));
tmp = find(tmp==0);
if ~isempty(tmp)
    start_time = tmp(end)+1;
    y = y(:,start_time:end-1);
end

inds_y = ones(size(y));
inds_y(find(y==0)) = 0;
inds_y = inds_y > 0;

[p N] = size(y);

x = [1:N]./N;

c = 100;
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
prior_params.hypers.a1 = 10;
prior_params.hypers.a2 = 10;

settings.L = 10;
settings.k = 20;
settings.Niter = 10000;
settings.saveEvery = 100;
settings.storeEvery = 10;
settings.saveMin = 1;
settings.saveDir = 'flu';
settings.trial = 1;
settings.init2truth = 0;
settings.sample_K_flag = 3;
settings.latent_mean = 1;
settings.inds_y = inds_y;

BNP_covreg_varinds(y,prior_params,settings,0);

