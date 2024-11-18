classdef Immodel
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        u
        v
        w
        F
        F_2D
        Fextrap
        Fextrap_2D
        uu
        vv
        ww
        xx
        yy
        F_im
        F_2D_im
        dFdu_im
        dFdv_im
        dFdw_im
        dFdu_2D_im
        dFdv_2D_im
        dFdu
        dFdv
        dFdw
        dFdu_2D
        dFdv_2D
    end
    
    methods
        function immodel = Immodel(template)
            %immodel should have fields
            %x (or u)
            %y (or v)
            %z (or w?) 
            %F(x,y,z) -- functional, not matrix form
            %F_2D(x,y)
            %Fextrap
            %
            %all x,y,z should be equally spaced
            
            
            if (isfield(template, 'x'))
                u = template.x;
            else
                u = template.u;
            end
            if (isfield(template, 'y'))
                v = template.y;
            else
                v = template.v;
            end
            if (isfield(template, 'z'))
                w = template.z;
            else
                w = template.w;
            end
            
            if (~isfield(template, 'F_2D'))
                template.F_2D = max(template.F, [], 3);
            end
            
            if ~isfield(template, 'Fextrap')
                template.Fextrap = min(template.F,[],[1 2 3]);
                template.Fextrap_2D = min(template.F_2D,[],[1 2 3]);
            end
            
            immodel.Fextrap = template.Fextrap;
            immodel.Fextrap_2D = template.Fextrap_2D;
            
            pdlen = 5;
            newu = interp1(1:length(u), u, (1-pdlen):(length(u)+pdlen), 'linear', 'extrap');
            newv = interp1(1:length(v), v, (1-pdlen):(length(v)+pdlen), 'linear', 'extrap');
            neww = interp1(1:length(w), w, (1-pdlen):(length(w)+pdlen), 'linear', 'extrap');
            immodel.u = newu;
            immodel.v = newv;
            immodel.w = neww;
            [immodel.uu, immodel.vv, immodel.ww] = ndgrid(newu, newv, neww);
            [immodel.xx, immodel.yy] = ndgrid(newu, newv);
            immodel.F_im = padarray(template.F, [pdlen pdlen pdlen], immodel.Fextrap);
            immodel.F_2D_im = padarray(template.F_2D, [pdlen pdlen], immodel.Fextrap_2D);


            deltau = median(diff(newu),'omitnan');
            deltav = median(diff(newv),'omitnan');
            deltaw = median(diff(neww),'omitnan');
            
  
             %simoncelli derivatives w/out 2nd order
             k = [.037659 .249153 .426375 .249153 .037659];
             d = [.109604 .276691 0 -.276691 -.109604];
  
             %most basic possible - works better with minimizer
             k = [0 0 1 0 0];
             d = [0 .5 0 -.5 0];
            
            
            immodel.dFdu_im = convn(convn(convn(padarray(immodel.F_im, [2 2 2], 'replicate'), d', 'valid'), k, 'valid'), reshape(k, 1,1,[]), 'valid') ./deltau;
            immodel.dFdv_im = convn(convn(convn(padarray(immodel.F_im, [2 2 2], 'replicate'), k', 'valid'), d, 'valid'), reshape(k, 1,1,[]), 'valid') ./deltav;
            immodel.dFdw_im = convn(convn(convn(padarray(immodel.F_im, [2 2 2], 'replicate'), k', 'valid'), k, 'valid'), reshape(d, 1,1,[]), 'valid') ./deltaw;
            
            immodel.dFdu_2D_im = convn(convn(padarray(immodel.F_2D_im, [2 2], 'replicate'), d', 'valid'), k, 'valid') ./deltau;
            immodel.dFdv_2D_im = convn(convn(padarray(immodel.F_2D_im, [2 2], 'replicate'), k', 'valid'), d, 'valid') ./deltav;
           
            
            
            fields = {'F', 'dFdu', 'dFdv', 'dFdw'};
           
            for j = 1:length(fields)
              immodel.(fields{j}) = griddedInterpolant(immodel.uu, immodel.vv, immodel.ww, immodel.([fields{j} '_im']), 'linear', 'nearest');
            end    
            fields = {'F_2D', 'dFdu_2D', 'dFdv_2D'};
          
            for j = 1:length(fields)
              immodel.(fields{j}) = griddedInterpolant(immodel.xx, immodel.yy, immodel.([fields{j} '_im']), 'linear', 'nearest');
            end         
            
            
        end
        
        function [g, gdelta] = gradientTest(immodel, pt, delta)
            if (nargin < 3 || isempty(delta))
                delta = 1e-3;
            end
            if (size(pt,1) == 2)
                g = [immodel.dFdu_2D(pt(1,:), pt(2,:));immodel.dFdv_2D(pt(1,:), pt(2,:))];
                gdelta = [immodel.F_2D(pt(1,:) + delta,pt(2,:)) - immodel.F_2D(pt(1,:), pt(2,:)); immodel.F_2D(pt(1,:), pt(2,:)+delta) - immodel.F_2D(pt(1,:), pt(2,:))]/delta;
            else
                 g = [immodel.dFdu(pt(1,:), pt(2,:),pt(3,:));immodel.dFdv(pt(1,:), pt(2,:),pt(3,:));immodel.dFdw(pt(1,:), pt(2,:),pt(3,:))];
                gdelta = [immodel.F(pt(1,:) + delta,pt(2,:),pt(3,:)) - immodel.F(pt(1,:), pt(2,:),pt(3,:)); immodel.F(pt(1,:) ,pt(2,:)+ delta,pt(3,:)) - immodel.F(pt(1,:), pt(2,:),pt(3,:));immodel.F(pt(1,:),pt(2,:),pt(3,:) + delta) - immodel.F(pt(1,:), pt(2,:),pt(3,:))]/delta;
            end
        end
        
    end
        
end
        

