function [tfit, mag] = energyMinTAlign (txknown, vknown, vunknown, tends, mindelta,  alpha_spring, tinitial, maginitial)
%function [tfit] = energyMinTAlign (txknown, vknown, vunknown, tends, mindelta,  alpha_spring, tinitial, maginitial)
%sigma for derivatives
%existsAndDefault('sigma', 0.01);


if (size(vknown,1) < size(vknown,2))
    vknown = vknown';
end
if (size(vunknown,1) < size(vunknown,2))
    vunknown = vunknown';
end


sigma = 0.01;
tx = txknown;
vs = vknown;
dvs = deriv(vs, sigma); %convolve with [0.5 0 -0.5] if sigma = 0.01
ddvs = deriv(dvs, sigma);
K = 1000;

if (size(dvs,1) < size(dvs,2))
    dvs = dvs';
end
if (size(ddvs,1) < size(ddvs,2))
    ddvs = ddvs';
end

if (size(vs,1) ~= length(txknown))
    error ('txknown and vknown must have same size');
end


% if (~existsAndDefault('tinitial', []))
%     tinitial = linspace(tends(1) + mindelta, tends(2) - mindelta, size(vunknown,1))';
%    % tinitial([1 end])
%     [xc,lags] = xcorr(interp1(txknown, vs(:,1),tinitial), vunknown(:,1), 10);
%     %plot(xc, lags); 
%     [~,I] = max(xc);
%     deltat = median(diff(tinitial));
%     tinitial = linspace(tends(1) + max(mindelta, lags(I)*deltat), tends(2) + min(-mindelta, lags(I)*deltat),  size(vunknown,1))';
%   %  tinitial([1 end])
% end

if(size(tinitial, 1) < size(tinitial,2))
    tinitial = tinitial';
end
if ~existsAndDefault('maginitial', 1)
    sk = sqrt(sum(vknown.^2,2));
    su = sqrt(sum(vunknown.^2,2));
    maginitial = mean(su, 'omitnan')/mean(sk(tx > tends(1) & tx < tends(2)), 'omitnan');
end
problem.objective = @totalEnergy;
problem.x0 = [maginitial; tinitial];
problem.solver = 'fminunc';
problem.options = optimoptions(problem.solver);
problem.options.GradObj = 'on';
problem.options.Hessian = 'on';
problem.options.Algorithm = 'trust-region';
problem.options.Display = 'none';
[testf, testg, testh] = problem.objective(problem.x0);
spd = spdiags(testh, [-1 1]);
inds = 1:(size(spd,1)-1);
plot (inds, spd(1:end-1,:), inds, spd(2:end, :))


%             testf
% bob = problem.objective(problem.x0 - 1e-6*testg./sqrt(sum(testg.^2))) -  testf
x = fminunc(problem); %, [], 2);
mag = x(1);
tfit = x(2:end);
%[~, aGrad, hessE1] = alignmentEnergy(t);
% [~, msg] = minSpacingEnergy(tfit);
% [~, sprg] = springEnergy(tfit);
% figure(1); plotyy(2:length(tfit), diff(tfit), 1:length(tfit), msg);
% figure(2); plotyy(2:length(tfit), diff(tfit), 1:length(tfit), sprg);


    function [E, gradE, hessE] = alignmentEnergy(x)
        t = x(2:end);
        m = x(1);
        dv = interp1(tx, vs, t, 'linear') - m*vunknown;
        E = 0.5 * sum(sum(dv.^2));
        gradE = [-sum(dot(dv,vunknown,2)); sum(dv.*interp1(tx,dvs,t,'linear'),2)];
        hessE = sum(dv.*interp1(tx,ddvs,t,'linear') + interp1(tx,dvs,t,'linear').^2,2);
        hessE = spdiags([0; hessE], 0, length(hessE)+1, length(hessE)+1);
        hess0 = [sum(sum(vunknown.^2)); -dot(vunknown,interp1(tx,dvs,t,'linear'),2)];
        hessE(1,:) = hess0;
        hessE(:,1) = hess0;
    end

    function [E, gradE, hessE] = minSpacingEnergy(x)
        t = x(2:end);
        tt = [tends(1); t; tends(end)];
        dt = exp(-K*(diff(tt) - mindelta)); %dt(1) = t(1)-t0 dt(2) = t(2) - t(1)
        dtf = dt(2:end);
        dtr = dt(1:(end-1)); %dtr(1) = t(1) - t0; dtf(1) = t(2)-t(1)
        E = sum(dt);
        gradE = [0;K*(dtf - dtr)];
        d2E = [0;K^2*(dtr - dtf)];
        d2Edn1 = [0;0;-K*dtr(1:end-1)];
        d2Edp1 = [0;-K*dtf(1:end)];
        hessE = spdiags([d2Edp1 d2E d2Edn1], -1:1,  length(d2E), length(d2E));
    end
    function [E, gradE, hessE] = springEnergy(x)
         t = x(2:end);
        tt = [tends(1); t; tends(end)];
        dt = diff(tt);
        E = 0.5*alpha_spring*sum(dt.^2);
        gradE = [0;alpha_spring*(dt(1:(end-1)) - dt(2:end))];
        n = length(t);
        e = alpha_spring*ones(n,1);
        hessE = spdiags([-[0;e] 2*[0;e] -[0;0;e(1:end-1)]], -1:1, n+1, n+1);
    end
    function [E, gradE, hessE] = totalEnergy(t)
        [E1, gradE1, hessE1] = alignmentEnergy(t);
        [E2, gradE2, hessE2] = minSpacingEnergy(t);
        [E3, gradE3, hessE3] = springEnergy(t);
        E =E1 + E2 + E3;
        gradE = gradE1 + gradE2 + gradE3;
        hessE = hessE1 + hessE2 + hessE3;
    end

            
end

