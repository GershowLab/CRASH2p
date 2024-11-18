function vr = parseHyperscopeDataFromBinary(srcdir)

fprintf('Parsing hyperscope data from %s...\n',srcdir);
tic;

vr = VolumeReconstructor(srcdir);
vr = vr.getStage();
tw_lowpass = 0.1;
tw_galvoDT = 5;
trange_galvoDT = vr.getStillTimes2(tw_lowpass,tw_galvoDT);
expectedGalvoDeltaT = -9.15e-5;
scanbeam = 2;
try
    vr = vr.calculateGalvoDeltaT(scanbeam,'a',trange_galvoDT);
    if abs((vr.galvoDeltaT-expectedGalvoDeltaT)/expectedGalvoDeltaT)>0.05
        % something wrong with galvoDeltaT calculation; use anticipated value instead
        vr.galvoDeltaT = expectedGalvoDeltaT;
    end
    vr = vr.calculatePhaseCorrection();
    vr = vr.pongFrames();
    vr = vr.applyScanOffset();
    vr = vr.addBehaviorVideo();
catch me
    disp(me.getReport());
    return;
end
% make more necessary calculations on vr.tsr
tsr = vr.tsr;
tsr = tsr.addBehaviorVideo();
try
    tsr = tsr.addSleapResult();
    tsr = tsr.findCounterMovements();
catch me
    warning('error loading SLEAP behavior labels; make sure SLEAP annotation is done and exported .h5 data is located in the behavior_video/ folder, then do tsr.addSleapResult() and tsr.findCounterMovements() later');
    disp(me.getReport);
end
vr.tsr = tsr;
vr.z_scale = vr.settings(end).analogScalingFactors(2).zTagMicronsPerTick*160E6/(2*pi*vr.settings(end).analogScalingFactors(2).tagLensFrequency);
% save parsed vr to disk
[~,timestamp,~] = fileparts(srcdir);
fpath = fullfile(srcdir,sprintf('%s_vr.mat',timestamp));
save(fpath,'vr','-v7.3');
fprintf('...done! Took %.1fs, data saved to %s\n',toc,fpath);

end