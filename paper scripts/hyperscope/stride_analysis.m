function [strideAlignedData,data] = stride_analysis(data, parameters)


if (nargin == 1 && isfield(data,'fwdaligned')) %fast redo if needed, mainly for debugging
    strideAlignedData = pca_metrics(data);
    return;
end

if (nargin == 1)
    expt = data;
    parameters = expt;
    srcdir = expt.srcdir;
    register = expt.register;
    sleap_path = expt.sleap_path;
    [~,timestamp,~] = fileparts(expt.srcdir);
    vr_name = [timestamp '_vr.mat'];
    data = load(fullfile(srcdir,vr_name));
    data.vr = data.vr.changeFilename(srcdir);
    data.tsr = data.vr.tsr.conditionData();
    data.tsr = data.tsr.addBehaviorVideo();
    try
        data.tsr = data.tsr.addSleapResult(sleap_path);
    catch
         disp('did not load sleap results');
         return;
    end
    data.tsr = data.tsr.findCounterMovements();
 
    
    template_name = [timestamp '_template.mat'];    
    s = load(fullfile(srcdir, register, template_name));
    data = mergeStructs(data,s);
    if (isfield(expt, 'projinds') && ~isempty(expt.projinds))
        expt.ladderinds = expt.projinds;
    end
    if (isfield(expt, 'ladderinds') && ~isempty(expt.ladderinds))
        zi = expt.ladderinds.zi;
        if (islogical(zi))
            z1 = find(zi,1);
            z2 = find(zi, 1, 'last' );
        else
            z1 = min(zi);
            z2 = max(zi);
        end
        fname = [timestamp sprintf('_z%dto%d_xytprojections.mat', z1,z2)];
    else
        fname = [timestamp '_xytprojections.mat'];
    end
   % try
        s = load(fullfile(srcdir, register, fname));
        data = mergeStructs(data,s);
    %catch
     %   disp ('xyt projections not loaded');
    %end
end



proj_sigma = [1 1 1];
gc2 = imblur(data.projections.gc, proj_sigma);
gd2 = imblur(data.projections.gd, proj_sigma);
rc2 = imblur(data.projections.rc, proj_sigma);
rd2 = imblur(data.projections.rd, proj_sigma);





f = {'fwdaligned','bckaligned'};
c = {'green', 'red','rat'};
tf = {'tfwd', 'tbck'};

counts = {gc2, rc2};
dwell = {gd2, rd2};

% counts = {gc2./sqrt(alpha),rc2./sqrt(alpha)};
% dwell = {gd2.*sqrt(alpha),rd2.*sqrt(alpha)};

if (~isfield(parameters,'pca_xrange') || isempty(parameters.pca_xrange))
    parameters.pca_xrange = [-Inf Inf];
end
if (~isfield(parameters,'pca_yrange') || isempty(parameters.pca_yrange))
    parameters.pca_yrange = [-Inf Inf];
end

if (~isfield(parameters,'rmax') || isempty(parameters.rmax))
    parameters.rmax = 2e7;
end
if (~isfield(parameters,'gmax') || isempty(parameters.gmax))
    parameters.gmax = 2e7;
end

lambdalim = [parameters.gmax, parameters.rmax];
strideAlignedData.template_mask= imopen(imclose(data.final_template.F_2D > median(data.final_template.F_2D,'all'),ones(7)),ones(5));
strideAlignedData.template_mask(1,:) = 0;
strideAlignedData.template_mask(end,:) = 0;
strideAlignedData.template_mask(:,1) = 0;
strideAlignedData.template_mask(:,end) = 0;


for i = 1:2
    
    strideAlignedData.(f{i}).xaxis = data.final_template.u;%(1:2:end);
    strideAlignedData.(f{i}).yaxis = data.final_template.v;%(1:2:end);
    strideAlignedData.(f{i}).phase_axis = deg2rad(7.5:15:360);    
    if (isempty(parameters.(tf{i})))
        strideAlignedData.(f{i}).nstride = 0;
        continue;
    end
    for j = 1:2
        [~,ca,~,ne] =  data.tsr.alignDataToMotionPhase(data.projections.tx, counts{j},strideAlignedData.(f{i}).phase_axis, parameters.(tf{i}));
        ca = ca;
        [~,da,~,ne] =  data.tsr.alignDataToMotionPhase(data.projections.tx, dwell{j},strideAlignedData.(f{i}).phase_axis, parameters.(tf{i}));
        da = da;
        ca(da <= 0) = 0;
        strideAlignedData.(f{i}).(c{j}) = min(squeeze(sum(ca,3,'omitnan'))./squeeze(sum(da,3, 'omitnan')),lambdalim(j));
        %strideAlignedData.(f{i}).(c{j})(~isfinite(strideAlignedData.(f{i}).(c{j}))) = 0;
        strideAlignedData.(f{i}).([c{j} '_counts']) = ca;
        strideAlignedData.(f{i}).([c{j} '_dwell']) = da;
        strideAlignedData.(f{i}).nstride = size(ca,3);
    end
    if (strideAlignedData.(f{i}).nstride == 0)
        continue;
    end
    npts = [ceil(diff(strideAlignedData.(f{i}).xaxis([1 end]))/20), 4];
    ic = IntensityCorrectorBSpline(sum(strideAlignedData.(f{i}).red_counts,[3 4],'omitnan')./sum(strideAlignedData.(f{i}).red_dwell,[3 4],'omitnan'), npts);
    alpha = 0 * strideAlignedData.(f{i}).red_counts;
    tic
    for m = 1:size(alpha,3)
        for n = 1:size(alpha,4)
            alpha(:,:,m,n) = ic.calculate_correction(strideAlignedData.(f{i}).red_counts(:,:,m,n),strideAlignedData.(f{i}).red_dwell(:,:,m,n));
        end
    end
    toc
    alpha = min(alpha,10);
    alpha = max(alpha,0.1);
    
    for j = 1:2
        ca = strideAlignedData.(f{i}).([c{j} '_counts'])./sqrt(alpha); 
        da = strideAlignedData.(f{i}).([c{j} '_dwell']).*sqrt(alpha);
        ca(da <= 0) = 0;
        da(da <= 0) = 0;
        strideAlignedData.(f{i}).(c{j}) = min(squeeze(sum(ca,3,'omitnan'))./squeeze(sum(da,3, 'omitnan')),lambdalim(j));
    end

    for j = 1:2
        rc = sum(strideAlignedData.(f{i}).red_counts,[3 4],'omitnan');
        strideAlignedData.(f{i}).red_sum = rc;
        rc(~isfinite(rc)) = 0;
        if (~any(rc > 0))
            continue;
        end
        rm = (rc > median(rc,'all'))&strideAlignedData.template_mask;
        strideAlignedData.(f{i}).red_mask = bwconvhull(rm) & imdilate(rm, ones(5)); %new expand within convex hull
        rd = sum(strideAlignedData.(f{i}).red_dwell,[3 4],'omitnan');
        strideAlignedData.(f{i}).red_mask = strideAlignedData.(f{i}).red_mask&(rd > percentile(rd(strideAlignedData.(f{i}).red_mask), .05)); %remove 5% of lowest sampled pixels
       
        nphi = length(strideAlignedData.(f{i}).phase_axis);
        X = reshape(strideAlignedData.(f{i}).(c{j}), [], nphi);
        mX =  mean(X,2,'omitnan');
        X = X - mX;
        Y = X(strideAlignedData.(f{i}).red_mask,:);
        Y = Y(all(isfinite(Y),2),:);
        if (isempty(Y) || size(Y,1) < 3)
            continue;
        end
        [coeff,~] = pca(Y, 'NumComponents',2,'Economy',true,'Centered',false);
        xc = xcorr(coeff(:,1), coeff(:,2));
        if (diff(xc(nphi + [0 1])) > 0)
            coeff = coeff(:,[2 1]);
        end
        score = X*coeff;
        strideAlignedData.(f{i}).([c{j} '_pca_z']) = reshape(score*[1;1i], size(rc));
        strideAlignedData.(f{i}).([c{j} '_pca_coeff']) = coeff;
       % strideAlignedData.(f{i}).([c{j} '_pca_reconstr']) = reshape(score*coeff')
    end
    xi = strideAlignedData.(f{i}).xaxis >= min(parameters.pca_xrange) & strideAlignedData.(f{i}).xaxis <= max(parameters.pca_xrange);
    yi = strideAlignedData.(f{i}).yaxis >= min(parameters.pca_yrange) & strideAlignedData.(f{i}).yaxis <= max(parameters.pca_yrange);
    strideAlignedData.(f{i}).pca_xinds = xi;
    strideAlignedData.(f{i}).pca_yinds = yi;
    strideAlignedData.(f{i}).pca_xaxis = strideAlignedData.(f{i}).xaxis(xi);
    strideAlignedData.(f{i}).pca_yaxis = strideAlignedData.(f{i}).yaxis(yi);
    for j = 1:2
        rc = sum(strideAlignedData.(f{i}).red_counts,[3 4],'omitnan');
        rc(~isfinite(rc)) = 0;
       
        if (~any(rc(xi,yi) > 0,'all'))
            continue;
        end
        nphi = length(strideAlignedData.(f{i}).phase_axis);
        X = reshape(strideAlignedData.(f{i}).(c{j})(xi,yi,:), [], nphi);
        X = X - mean(X,2);
        Y = X;
        Y = Y(all(isfinite(Y),2),:);
        if (isempty(Y) || size(Y,1) < 3)
            continue;
        end
        [coeff,~] = pca(Y, 'NumComponents',2,'Economy',true,'Centered',false);
        xc = xcorr(coeff(:,1), coeff(:,2));
        if (diff(xc(nphi + [0 1])) > 0)
            coeff = coeff(:,[2 1]);
        end
        score = X*coeff;
        strideAlignedData.(f{i}).([c{j} '_pca_restricted_z']) = reshape(score*[1;1i], size(rc(xi,yi)));
        strideAlignedData.(f{i}).([c{j} '_pca_restricted_coeff']) = coeff;
    end
end
strideAlignedData = pca_metrics(strideAlignedData);
strideAlignedData.nfwd = strideAlignedData.fwdaligned.nstride;
strideAlignedData.nbck = strideAlignedData.bckaligned.nstride;


c = data.projections.gc;
d = data.projections.gd;
nmedian = ceil(10/median(diff(data.projections.tx)));
lambdamed = medfilt1(c./d,nmedian,[],3,'omitnan');
c = c-d.*lambdamed;
[~,ca] = data.tsr.alignDataToMotionPhase(data.projections.tx, c,strideAlignedData.fwdaligned.phase_axis, parameters.tfwd);
[~,da] = data.tsr.alignDataToMotionPhase(data.projections.tx, d,strideAlignedData.fwdaligned.phase_axis, parameters.tfwd);
%strideAlignedData.fwdaligned.dleta_green_counts = ca;
strideAlignedData.fwdaligned.delta_green = squeeze(sum(ca,3,'omitnan'))./squeeze(sum(da,3, 'omitnan'));

end

function strideAlignedData = pca_metrics(strideAlignedData)

    f = {'fwdaligned','bckaligned'};
    c = {'green', 'red'};
    for i = 1:2
        if (strideAlignedData.(f{i}).nstride <= 0)
            continue;
        end
        try
            [xx,yy] = ndgrid(strideAlignedData.(f{i}).xaxis, strideAlignedData.(f{i}).yaxis);
        catch
            continue;
        end
        x = xx(strideAlignedData.(f{i}).red_mask);
        y = yy(strideAlignedData.(f{i}).red_mask);
        nphi = length(strideAlignedData.(f{i}).phase_axis);
        for j = 1:2
            try
                z = strideAlignedData.(f{i}).([c{j} '_pca_z'])(strideAlignedData.(f{i}).red_mask);
            catch
                continue;
            end
%             [xa,mz_x] = meanyvsx(x,z,strideAlignedData.(f{i}).xaxis);
%             valid = isfinite(mz_x) & isfinite(xa);
%             xa = xa(valid);
%             mz_x = mz_x(valid);
%             phi = unwrap(angle(mz_x));
%            % mp = mean(phi);
%             %phi = phi - mp;
%            
%             p = polyfit(xa,phi,1);
            
            [k,phi0,~] = line_fit_phase_z_vnc(x, y, z);

            strideAlignedData.(f{i}).([c{j} '_pca_kvector']) =k;
            coeff = strideAlignedData.(f{i}).([c{j} '_pca_coeff']);
            u = mean(strideAlignedData.(f{i}).(c{j}),3,'omitnan');

            strideAlignedData.(f{i}).([c{j} '_pca_wave_z']) = reshape(abs(strideAlignedData.(f{i}).([c{j} '_pca_z'])).*exp(1i*(k(1)*xx + k(2)*yy +phi0)), size(u));
            zz = strideAlignedData.(f{i}).([c{j} '_pca_wave_z']); 
            strideAlignedData.(f{i}).([c{j} '_pca_wave_reconstr']) = real(zz).*reshape(coeff(:,1),1,1,[]) + imag(zz).*reshape(coeff(:,2),1,1,[]) + u;
            zz = strideAlignedData.(f{i}).([c{j} '_pca_z']);
            strideAlignedData.(f{i}).([c{j} '_pca_reconstr']) =  real(zz).*reshape(coeff(:,1),1,1,[]) + imag(zz).*reshape(coeff(:,2),1,1,[]) + u;

            rm = reshape(strideAlignedData.(f{i}).red_mask,[],1);
            X = reshape(strideAlignedData.(f{i}).(c{j}), [], nphi);
            X = X(rm,:);
            X_PCA = reshape(strideAlignedData.(f{i}).([c{j} '_pca_reconstr']),[],nphi);
            X_PCA = X_PCA(rm,:);
            X_PCA_W = reshape(strideAlignedData.(f{i}).([c{j} '_pca_wave_reconstr']),[],nphi);
            X_PCA_W = X_PCA_W(rm,:);
            mX = mean(X,2,'omitnan');

            strideAlignedData.(f{i}).([c{j} '_var_expl_pca']) = 1 - var(X-X_PCA,0,"all")./var(X-mX,0,'all');
            strideAlignedData.(f{i}).([c{j} '_var_expl_pca_wave']) = 1 -  var(X-X_PCA_W,0,'all')./var(X-mX,0,'all');
            
        end
    end
end


function a = mergeStructs(a,b)
    fn = fieldnames(b);
    for i = 1:length(fn)
        a.(fn{i}) = b.(fn{i});
    end
end
