function [theta,w,thetas,ws, zp, zps] = locationEstimateKalman (tp, xp, D, tx, measErrPerPhoton, g)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

   
    ndim = size(xp, 1);
    existsAndDefault('g', eye(ndim));
    w = zeros([2*ndim 2*ndim length(tx)]);
    wkkm1 = w;
    F = w;
    theta = zeros([2*ndim length(tx)]);
    thetakkm1 = theta;


    initialLoc = mean(xp(:,tp < tx(3)), 2);
    
    lastT = 2*tx(1)- tx(2);
    tt = [initialLoc; 0*initialLoc];
    ww = diag(1e6*ones([size(w,1),1]));
    if (any(size(D) == 1))
        D = diag(D);
    end
    ts1 = tic;
    for j = 1:length(tx)
        p1 = bsearch(tp, lastT);
        p2 = bsearch(tp, tx(j));
        deltaT = tx(j) - lastT;
        F(:,:,j) =  eye(2*ndim) + deltaT*diag(ones([1 ndim]), ndim);
        Q = deltaT*F(:,:,j)*D*F(:,:,j)';
        deltaTphoton = tp(p1:p2)-lastT;
        [thetakkm1(:,j), theta(:,j), wkkm1(:,:,j), w(:,:,j)] = singleStep (tt, ww, xp(:,p1:p2), deltaTphoton, F(:,:,j), Q, measErrPerPhoton, g);
        lastT = tx(j);
        ww = w(:,:,j);
        tt = theta(:,j);
        if (toc(ts1) > 10)
            plot(tx(1:j), theta(1:3,1:j)); xlim([tx(1) tx(end)]); pause(0.0001);
            ts1 = tic;
        end
    end
    
    ws = w;
    thetas = theta;
    
    for j = (length(tx)-1):-1:1
        C = w(:,:,j)*F(:,:,j+1)'/wkkm1(:,:,j+1);
        thetas(:,j) = theta(:,j) + C*(thetas(:,j+1) - thetakkm1(:,j+1));
        ws(:,:,j) = w(:,:,j) + C * (ws(:,:,j+1) - wkkm1(:,:,j+1)) * C';
    end
    zp = xp - interp1(tx, theta(1:ndim,:)', tp)';
    zps =  xp - interp1(tx, thetas(1:ndim,:)', tp)';
    
%     
%     for j = (length(tx)-1):-1:1
%         H = w(:,:,j)/wkkm1(:,:,j+1);
%         thetas(:,j) = theta(:,j) + H*(thetas(:,j+1) - thetakkm1(:,j+1));
%         ws(:,:,j) = w(:,:,j) + H * (ws(:,:,j+1) - wkkm1(:,:,j+1)) * H';
%     end
    
    
end

%ND
function [thetakkm1, thetakk, wkkm1, wkk] = singleStep (theta, w, xp, deltaTphoton, F, Q, r, g)

    ndim = size(xp,1);
    thetakkm1 = F*theta;
    wkkm1 = F*w*F' + Q;

    
    dxp = xp-repmat(theta(1:ndim),[1,size(xp,2)])-repmat(deltaTphoton,[ndim,1]).*theta(ndim + (1:ndim));
    validRangeMultiplier = 6;
    
    %numphotons = length(dxp);
    if (any(size(r) == 1))
        if (isrow(r))
            r = r';
        end
        rtest = repmat(validRangeMultiplier*r, [1 size(dxp,2)]);
    else
         rtest = repmat(validRangeMultiplier*diag(r), [1 size(dxp,2)]);
    end
    
    pvalid = all(dxp.^2 < rtest, 1);
    
    %if less than 10% of photons are valid, only restrict by Z
    %if that fails, then do not restrict
    if (mean(pvalid) < 0.1)
        pvalid = dxp(end,:).^2 < rtest(end,:);
    end
    if (mean(pvalid) < 0.1)
        pvalid = true(size(deltaTphoton));
    end
    numphotons = nnz(pvalid);
    if ~(any(pvalid))
        thetakk = thetakkm1;
        wkk = wkkm1;
        return;
    end
    
    dxp = dxp(:,pvalid);

    H = [diag(ones([1 ndim])) diag(zeros([1 ndim]))];

    if (any(size(r) == 1))
        r = diag(r);
    end
    if (any(size(g) == 1))
        g = diag(g);
    end
    
    y = sum(g*dxp,2)/numphotons;
    S = H*wkkm1*H' + r/numphotons;
    K = (wkkm1*H')/S;
    
    thetakk = thetakkm1 + K*y;
    wkk = wkkm1 - K*H*wkkm1;
    if (any(~isfinite(thetakk)) || any(~isfinite(wkk(:))))
        wkk = wkkm1;
        thetakk = thetakkm1;
        warning ('something foobar');
    end
        
end


% 1D
% function [thetakkm1, thetakk, wkkm1, wkk] = singleStep (theta, w, xp, deltaTphoton, F, Q, r, g)
% 
%     thetakkm1 = F*theta;
%     wkkm1 = F*w*F' + Q;
% 
%     dxp = xp-theta(1)-deltaTphoton*theta(2);
%     validRangeMultiplier = 10;
%     
%     %numphotons = length(dxp);
%     
%     pvalid = dxp.^2 < validRangeMultiplier*r;
%     numphotons = nnz(pvalid);
%     if ~(any(pvalid))
%         thetakk = thetakkm1;
%         wkk = wkkm1;
%         return;
%     end
%     
%     dxp = dxp(pvalid);
% 
%     d = numphotons*wkkm1(1,1) + r;
%     thetakk = thetakkm1 + g*[wkkm1(1,1);wkkm1(1,2)]*sum(dxp)/d;
%     wkk = wkkm1*r/d;
%     wkk(2,2) = wkkm1(2,2)-numphotons*(wkkm1(2,1).^2)/d;
%     
% end
