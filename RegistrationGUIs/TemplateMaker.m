classdef TemplateMaker < handle
    % helper class for TemplateCreator app
    %   contains input data and space for return values
    %   input: vra post alignment, must have functioning vr (volume reconstructor) to reload
    %       raw image data
    %   output: roi (fields x,y,z): the limits in microns of the 3D rectangle roi from
    %       which to make the template and proceed with non-rigid registration
    %           xybin: x&y dimension in microns of voxel for binning (typical value, 1)           
    %           zbin: z dimension in microns of voxel for binning (typical value, 2)
    %           smooth_s: setting for 1D spline smoother and interpolator
    %           selected_frames: which frames should be used to construct
    %           the template
    %           fa: fast assembler that can be used to make images
    %               compatible with template
    %           template: standard 3D template volume for use in macros and
    %                     scripts - constructed using fa, selected frames,
    %                     and smooth_s
    %   finished: a boolean flag the TemplateCreator app sets to true when
    %       output fields are filled in and we are ready to proceed
    %
    %  typical use:
    %           tm = TemplateMaker(vra)
    %           TemplateCreator(tm);
    %           waitfor(tm, 'finished', true);
    %           %return values are now stored in tm
    properties
        %inputs
        vra
        
        %outputs
        roi
        xybin
        zbin
        smooth_s
        method = 'Best'
        selected_frames

        %derived outputs (for convenience, determined by above)
        template
        fa

        %analysis mask outputs
        analysis_mask

        frame_buffer

        finished = false;
    end

    methods


        function tm = TemplateMaker(vra)
            tm.vra = vra;
        end

        function setBinning(tm, roi, xybin, zbin)
            try
                if (xybin == tm.xybin && zbin == tm.zbin && all(roi.x == tm.roi.x & roi.y == tm.roi.y & roi.z == tm.roi.z))
                    return;
                end
            catch
            end

            midx = mean(roi.x);
            midy = mean(roi.y);
            midz = mean(roi.z);
            nx = ceil(diff(roi.x)/xybin);
            ny = ceil(diff(roi.y)/xybin);
            nz = ceil(diff(roi.z)/zbin);

            xe = midx + ((-nx):2:(nx))*xybin/2;
            ye = midy + ((-ny):2:(ny))*xybin/2;
            ze = midz + ((-nz):2:(nz))*zbin/2;
           
% 
%             xc = (roi.x(1)+xybin/2):xybin:(roi.x(2)+xybin/2);
%             yc = (roi.y(1)+xybin/2):xybin:(roi.y(2)+xybin/2);
%             zc = (roi.z(1)+zbin/2):zbin:(roi.z(2)+zbin/2);
            
            tm.fa = FastAssembler(xe,ye,ze,[],false,tm.vra.z_scale);
            tm.roi = roi;
            tm.xybin = xybin;
            tm.zbin = zbin;
            
            fb.ind = 0;
            fb.counts = NaN(tm.fa.getImSize);
            fb.dwell = fb.counts;

            tm.frame_buffer = repmat(fb, [0 1]);
        end

        function [counts, dwell] = binFrames(tm, inds)
            counts = zeros(tm.fa.getImSize);
            dwell = zeros(tm.fa.getImSize);
            for j = 1:length(inds)
                fbi = find([tm.frame_buffer.ind] == inds(j), 1, 'first');
                if (isempty(fbi))
                    fb.ind = inds(j);
                    [fb.counts, fb.dwell] = tm.vra.binFrame(inds(j),tm.fa);
                    tm.frame_buffer = [tm.frame_buffer fb]; 
                else
                    fb = tm.frame_buffer(fbi);
                end
                counts = counts+fb.counts;
                dwell = dwell+fb.dwell;
            end
        end
% 
%         function im = makeImGauss(tm,inds,sigma)
% 
%             sigma = sigma./[tm.fa.dx tm.fa.dy tm.fa.dz];
%             [counts,dwell] = tm.binFrames(inds);
%             im = imblur(counts,sigma)./imblur(dwell,sigma);
%             mmat = imclose(dwell > 0, ones([9 9 9]));
%             im = mmat.*im;
%             im(im < tm.vra.rmin) = tm.vra.rmin;
% 
% 
%         end


        function im = makeImSpline(tm,inds,s)
            tm.smooth_s = s;
            tm.selected_frames = inds;
            [counts,dwell] = tm.binFrames(inds);
            im = l1spline(counts./dwell, s, 1, 100, 1, 1e-2);
            mmat = imclose(dwell > 0, ones([9 9 9]));
            im = mmat.*im;
            im(im < tm.vra.rmin) = tm.vra.rmin;
        end

        function str = toText(tm)
            str{1} = sprintf('roi x: %.2f\t%.2f\nroi y: %.2f\t%.2f\nroi z: %.2f\t%.2f\n', tm.roi.x(1), tm.roi.x(2), tm.roi.y(1), tm.roi.y(2), tm.roi.z(1), tm.roi.z(2));
            str{2} = sprintf('dx = dy = %.2f \t dz = %.2f\n', tm.xybin, tm.zbin);            
            str{3} = sprintf('smooth = %.4g\n', tm.smooth_s);
            str{4} = ['frames: ' sprintf('%d, ', tm.selected_frames)];
            str = [str{:}];
        end

        function makeAnalysisMask(tm, blur_sigma, threshold,dilate_px, close_px)
            tm.analysis_mask.sigma = blur_sigma;
            tm.analysis_mask.threshold = threshold;
            tm.analysis_mask.dilate = dilate_px;
            tm.analysis_mask.close = close_px;

            tm.analysis_mask.blurred_im = imblur(tm.template.F, blur_sigma./([tm.xybin tm.xybin tm.zbin]));
            mask = tm.analysis_mask.blurred_im > threshold;
            tm.analysis_mask.thresholded_mask = mask;
            if (dilate_px > 0)
                mask = imdilate(mask, strel('sphere',ceil(dilate_px)));
            end
            if (close_px > 0)                  
                mask = imclose(mask, strel('sphere', ceil(close_px)));
            end
            tm.analysis_mask.mask = mask;


        end

        function templateSettings = getSettingsStruct(tm)
            fn = {'roi', 'xybin', 'zbin', 'smooth_s', 'selected_frames', 'fa', 'template', 'analysis_mask'};
            for j = 1:length(fn)
                templateSettings.(fn{j}) = tm.(fn{j});
            end
        end


    end
end