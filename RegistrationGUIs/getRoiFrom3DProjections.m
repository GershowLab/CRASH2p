function roi = getRoiFrom3DProjections(xedges, yedges, zedges, volume, label)
%function roi = getRoiFrom3DProjections(xedges, yedges, zedges, volume, label)
%   wrapper for RoiSelectorInput.getRoi
%   displays an interactive gui and lets user select ROI in 3 projections
    existsAndDefault('label',{});
    roi = RoiSelectorInput.getRoi(xedges, yedges, zedges, volume, label);
end