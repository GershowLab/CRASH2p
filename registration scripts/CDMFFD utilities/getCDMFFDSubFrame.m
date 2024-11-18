function [counts,dwell] = getCDMFFDSubFrame(vra,fa,frameind,cdmffd_shifts,altcolor)

% get rigid transforms
time_edges = vra.subframe_time_edges(frameind,:);
nsubframes = length(time_edges)-1;
pvec = squeeze(vra.pvec_subframe(frameind,:,:));
if any(vra.additional_rotation)
    rotation_angle_rad = deg2rad(vra.additional_rotation([3 2 1])); % match pvec order
    rotation_angle_rad(3) = -rotation_angle_rad(3); % match image display axes direction
    pvec_add = [rotation_angle_rad 0 0 0];
    at_add = RigidAligner3D.rigidTransform(pvec_add);
    at_add(4,:) = [0 0 0 1];
end

% apply rigid and non-rigid transforms to photon and dwell locations
counts = zeros([fa.getImSize nsubframes]);
dwell = counts;
for k=1:nsubframes
    te = time_edges(k+[0 1]);
    if any(vra.additional_rotation)
        at = RigidAligner3D.rigidTransform(pvec(k,:));
        at(4,:) = [0 0 0 1];
        att = at_add*at; % order of operation: right to left
        att = att(1:3,:);
    else
        att = RigidAligner3D.rigidTransform(pvec(k,:));
    end
    % get raw photon and dwell locations
    [ppos,pphase] = vra.getPhotonLocations(te,altcolor);
    [dpos,dphase,~,dtau] = vra.getDwellLocationsSubsampled(te,[],altcolor);
    % convert to locvec (throw out old stage z and add phase-calculated z)
    plocvec = ppos;
    plocvec(3,:) = vra.z_scale*cos(pphase);
    dlocvec = dpos;
    dlocvec(3,:) = vra.z_scale*cos(dphase);
    % apply rigid transform to locvec
    plocvec = att*[plocvec;ones(1,size(plocvec,2))];
    dlocvec = att*[dlocvec;ones(1,size(dlocvec,2))];
    % apply non-rigid transform to locvec then bin image
    [ct,dw] = binWithCDMFFDShifts(plocvec,pphase,dlocvec,dtau,fa,cdmffd_shifts);
    counts(:,:,:,k) = ct;
    dwell(:,:,:,k) = dw;
end

end