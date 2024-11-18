function [templateSettings, tm] = setupTemplate(vra)
% function [templateSettings, tm] = setupTemplate(vra)
% vra < VRAligner post alignment with working vr
% or 
%  [templateSettings, tm] = setupTemplate(tm)
% tm < TemplateMaker class
% templateSettigns - has settings to make template in structure that's
% easily saved to disk
%
% wrapper for template maker and template creator
    if (isa(vra, 'TemplateMaker'))
        tm = vra;
    else
        tm = TemplateMaker(vra);
    end
    TemplateCreator(tm);
    waitfor(tm, 'finished', true)
    templateSettings = tm.getSettingsStruct();
end