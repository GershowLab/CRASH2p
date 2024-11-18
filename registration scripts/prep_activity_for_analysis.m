function activity = prep_activity_for_analysis(activity, trange, lpsigma)
    existsAndDefault('trange', [min(activity.tx) max(activity.tx)]);
    existsAndDefault('lpsigma',1);
   

    dwellcorrect = true;

    if (~isfield(activity.voisets(1), 'gdfine_ic'))
        warning('rerun activity_in_vois to capture dwell intensity correction -- using counts for now');
        dwellcorrect = false;
    end

    try
        %changed 6/19/2023 by MHG to use intensity corrected values
        for j = 1:length(activity.voisets)
            if (dwellcorrect)
                activity.voisets(j).ratio_fine =  activity.voisets(j).gcfine.* activity.voisets(j).rdfine_ic./( activity.voisets(j).gdfine_ic.* activity.voisets(j).rcfine);
            else
                activity.voisets(j).ratio_fine =  activity.voisets(j).gcfine_ic.* activity.voisets(j).rdfine./( activity.voisets(j).gdfine.* activity.voisets(j).rcfine_ic);
            end
        end
        for j = 1:length(activity.voisets)
            if (dwellcorrect)
                activity.voisets(j).lpr = lowpass1D((activity.voisets(j).gcfine.* activity.voisets(j).rdfine_ic)',lpsigma)'./lowpass1D((activity.voisets(j).gdfine_ic.* activity.voisets(j).rcfine)', lpsigma)';%lowpass1D(activity.voisets(j).ratio_fine',lpsigma)';
            else
                activity.voisets(j).lpr = lowpass1D((activity.voisets(j).gcfine_ic.* activity.voisets(j).rdfine)',lpsigma)'./lowpass1D((activity.voisets(j).gdfine.* activity.voisets(j).rcfine_ic)', lpsigma)';
            end
        end
    catch
        warning('Rerun activity_in_vois to capture dwell intensity correction -- using uncorrected intensities for now')
        for j = 1:length(activity.voisets)
            activity.voisets(j).ratio_fine =  activity.voisets(j).gcfine.* activity.voisets(j).rdfine./( activity.voisets(j).gdfine.* activity.voisets(j).rcfine);           
        end
        for j = 1:length(activity.voisets)
            activity.voisets(j).lpr = lowpass1D((activity.voisets(j).gcfine.* activity.voisets(j).rdfine)',lpsigma)'./lowpass1D((activity.voisets(j).gdfine.* activity.voisets(j).rcfine)', lpsigma)';%lowpass1D(activity.voisets(j).ratio_fine',lpsigma)';
        end
    end

    valid = activity.txfine >= min(trange) & activity.txfine <=max(trange);

    for j = 1:length(activity.voisets)
        %find baseline ratio
        activity.voisets(j).r0 = zeros([size(activity.voisets(j).lpr,1) 1]);
        for k = 1:size(activity.voisets(j).lpr,1)
            activity.voisets(j).r0(k) = meanAndStdAssymOutlier(activity.voisets(j).lpr(k,valid));
        end
    end