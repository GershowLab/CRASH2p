function fitstruct = photonSPPFtworates(et, npr, npg, d_loglambda, presmooth)
%function [lambda, lambda_s] = photonSPPFsinglerate(et, np,  d_loglambda)
% et is observation time, npr is number of red photons, npg is number of
% green photons
% et is either 1 x N or 2 x N
% if et is 1 x N, photon counting time (dt_obs) = diff(et)
% otherwise dt_obs = et(2,:) & et = et(1,:)
% uncertainty always increases with diff(et);

% fits a model lambda_red = exp(theta1); lambda_green = exp(theta1 +
% theta2)
% thus exp(theta2) is the ratio of green to red, a proxy for neural
% activity if green channel is gcamp signal and red channel is reference
%
% d_loglambda is diffusion constant in log(lambda), that is
% p(log(lambda(t+dt)) - log(lambda(t)) ~ N(0, d_loglambda dt)
% d_loglambda = 1 is a pretty reasonable starting point
% diffusion constant for the ratio is set at 0.1 * diffusion constant of
% the rate
%
% presmooth, if an integer greater than 0 mixes each observation with the
% (n = presmooth) points before and after it before running through the
% filter; this avoids overfitting to anti-correlated errors (e.g. due to
% bad photon binning)
%

%get initial guess to be a reasonable number
buffelems = 100;
existsAndDefault('presmooth', 0);
if (length(et) > buffelems)
    fs1 = photonSPPFtworates(et(:,1:buffelems), npr(1:buffelems), npg(1:buffelems), d_loglambda, presmooth);
    precalc = true;
else
    precalc = false;
end
ii = find(npr > 0 & npg > 0, 1, 'first');
et = et(:,ii:end);
npr = npr(ii:end);
npg = npg(ii:end);

if (size(et,1) > 1)
    dt_obs = et(2,:);
    et = et(1,:);
    dt = diff(et);
else
    dt = diff(et);
    dt_obs = [dt(1) dt];
end
dt = [dt(1) dt];

presmooth = round(presmooth);
if (presmooth >= 0)
    %resample each to have the average number of photons and observation
    %time from its presmooth nearest neighbors
    ck = ones([1 2*presmooth+1])/(2*presmooth+1);
    gw = npg.*dt_obs;
    rw = npr.*dt_obs;
    dt_obs = conv(dt_obs, ck, 'same');
    npg = conv(gw, ck, 'same')./dt_obs;
    npr = conv(rw, ck, 'same')./dt_obs;
end
w = zeros(2,2,length(et));
wjjm1 = w;
theta = zeros(2,length(et));
innov = theta;

if (precalc)
    theta(:,1) = fs1.theta(:,end);
    w(:,:,1) = fs1.w(:,:,end);
else
    
    theta(1,1) = log(max(1,npr(1))/dt_obs(1));
    theta(2,1) = log(max(1,npg(1))/max(1,npr(1)));
    w(:,:,1) = 1/(max(1,npg(1)));
    w(1,1,1) = 1/(max(1,npr(1) + npg(1)));
end
wjjm1(:,:,1) = w(:,:,1);

for j = 2:length(theta)
    if (numel(d_loglambda) == 1)
        q = d_loglambda * dt(j)* [sqrt(2) 0; 0 1/sqrt(2)];% * [1 0; 0 0.1];
    else
        q = d_loglambda * dt(j);
    end
    wjjm1(:,:,j) = w(:,:,j-1) + q;
    lgdt = exp(sum(theta(:,j-1))) * dt_obs(j);
    lrdt = exp(theta(1,j-1))*dt_obs(j);
    dwi = [lgdt+lrdt, lgdt;lgdt, lgdt];
    w(:,:,j) = inv(inv(wjjm1(:,:,j)) + dwi);
    if (any(any (~isfinite(w(:,:,j)))) || det(w(:,:,j)) == 0)
        warning ('w matrix has an issue in sppf');
        break;
    end
    innov(:,j) = [npr(j) - lrdt; npg(j) - lgdt];
    theta(:,j) = theta(:,j-1) + w(:,:,j) * [npr(j) + npg(j) - lgdt - lrdt; npg(j) - lgdt];
end

lambda_red = exp(theta(1,:));
lambda_green = exp(theta(1,:) + theta(2,:));
slambda_red = lambda_red.*sqrt(squeeze(w(1,1,:))).';
slambda_green = lambda_green.*sqrt(squeeze(sum(sum(w)))).';
ratio = exp(theta(2,:));
sratio = ratio.*sqrt(squeeze(w(2,2,:))).';
%implements smoothing step from Koyama, S., Eden, U.T., Brown, E.N., and Kass, R.E. (2009). Bayesian decoding of neural spike trains. Ann Inst Stat Math 62, 37.
%eqn 17-19
%assuming F = 1
theta_s = theta;
ws = w;

for j = (length(w)-1):-1:1
     if (any(any (~isfinite(w(:,:,j)))) || det(w(:,:,j)) == 0)
        warning ('w matrix has an issue in sppf');
        break;
     end
    h = w(:,:,j)/wjjm1(:,:,j+1);
    theta_s(:,j) = theta(:,j) + (h*(theta_s(:,j) - theta(:,j)));
    ws(:,:,j) = ws(:,:,j) + h*(ws(:,:,j+1)-wjjm1(:,:,j+1))*h';
end

lambda_red_s = exp(theta_s(1,:));
lambda_green_s = exp(theta_s(1,:) + theta_s(2,:));
slambda_red_s = lambda_red_s.*sqrt(squeeze(ws(1,1,:))).';
slambda_green_s = lambda_green_s.*sqrt(squeeze(sum(sum(ws)))).';
ratio_s = exp(theta_s(2,:));
sratio_s = ratio_s.*sqrt(squeeze(ws(2,2,:))).';


fn = {'et', 'theta', 'theta_s', 'w', 'ws', 'lambda_red', 'lambda_green', 'slambda_red', 'innov', ...
    'slambda_green', 'ratio', 'sratio','lambda_red_s', 'lambda_green_s', 'slambda_red_s', 'slambda_green_s', 'ratio_s', 'sratio_s'};
for j = 1:length(fn);
    fitstruct.(fn{j}) = eval(fn{j});
end