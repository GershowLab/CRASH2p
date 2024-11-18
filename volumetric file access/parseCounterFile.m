function [counterData, explanations] = parseCounterFile(fname)
% function [counterData, explanations] = parseCounterFile(fname)
% parses [date code]_str_ctrs (streaming counters) file and does minor
% error checking
%
% fname < string : path to file
% counterData < struct : parsed file; rollover includes only times when
%   160 MHz counter increased (other lines indicate times when number of
%   elements received increased)
% explanations < cell array of strings : explanation of what the fields
% mean
    data = importdata(fname);
    ch = matlab.lang.makeValidName(data.colheaders);
    fid = fopen(fname, 'rt');
    fgetl(fid);    
    raw.mat = fscanf(fid, '%lu', size(data.data'))';
    fclose(fid);
    for j = 1:length(ch)
        raw.(ch{j}) = raw.mat(:,j);
    end
    valid = raw.beam1_photon_u16s_in_to_pmt_volumetric > 0 | raw.beam1_tag_triggers_into_pmt_volumetric > 0 | raw.total_u16s_sent_to_encode_timers_and_aod_signals > 0 | ...
        raw.beam2_photon_u16s_in_to_pmt_volumetric > 0 | raw.beam2_tag_triggers_into_pmt_volumetric > 0 | raw.bytes_into_send_volumetric_as_u16 > 0 | raw.elements_received_from_fpga > 0;
    for j = 1:length(ch)
       raw.(ch{j}) = raw.(ch{j})(valid,:)';
    end
    raw.mat = raw.mat(valid,:);
    
    counterData.raw = raw;
    
    counterData.allCounterValues = raw.pmt_loop_counter___160_MHz/2^28;
    counterData.allCountersOk = ~(any(any(bitand(raw.pmt_loop_counter___160_MHz, uint64(2^28-1))))); %if this is true, this means one of the counters is not a multiple of 2^28; only supposed to update on rollovers
    counterData.rolloverMat = raw.mat(diff(counterData.allCounterValues) ~= 0,:); %values when the counter incremented, as opposed to when counters incremented when data was received
    for j = 1:length(ch)
       counterData.rollover.raw.(ch{j}) = counterData.rolloverMat(:,j)';
    end
    %counter number and data stream
    counterData.rollover.cnum = counterData.rollover.raw.pmt_loop_counter___160_MHz/2^28;
    counterData.rollover.bytesSentToComputer = diff(counterData.rollover.raw.dma_loop_bytes_sent_to_dma_fifo([1:end end])); %may not be the same number of bytes as in the file b/c of buffering
    counterData.rollover.bytesInCounter = diff(counterData.rollover.raw.bytes_into_send_volumetric_as_u16([1:end end]));
    
    %photons detected by the PMT, independent measure from volumetric loop
    counterData.rollover.photonsDetectedb1c0 = diff(counterData.rollover.raw.b1ch0_num_pmt_signals_integrated([1:end end]));
    counterData.rollover.photonsDetectedb1c1 = diff(counterData.rollover.raw.b1ch1_num_pmt_signals_integrated([1:end end]));
    counterData.rollover.photonsDetectedb2c0 = diff(counterData.rollover.raw.b2ch0_num_pmt_signals_integrated([1:end end]));
    counterData.rollover.photonsDetectedb2c1 = diff(counterData.rollover.raw.b2ch1_num_pmt_signals_integrated([1:end end]));
    
    %photons processed by volumetric loop
    counterData.rollover.photons_b1 = diff(counterData.rollover.raw.beam1_photon_u16s_in_to_pmt_volumetric([1:end end]));
    counterData.rollover.photons_b2 = diff(counterData.rollover.raw.beam2_photon_u16s_in_to_pmt_volumetric([1:end end]));

    %photons processed by volumetric loop
    counterData.rollover.tags_b1 = diff(counterData.rollover.raw.beam1_tag_triggers_into_pmt_volumetric([1:end end]));
    counterData.rollover.tags_b2 = diff(counterData.rollover.raw.beam2_tag_triggers_into_pmt_volumetric([1:end end]));

    %photons detected by pmt but not sent to volumetric loop
    %this can be OK if microscope was not sending data, e.g. if it's idle
    %if idle, num unprocessed photons should be about 1.6777 * dark rate
    counterData.rollover.photonsUnprocessed_b1 =  counterData.rollover.photonsDetectedb1c0 +  counterData.rollover.photonsDetectedb1c1 - counterData.rollover.photons_b1;
    counterData.rollover.photonsUnprocessed_b2 =  counterData.rollover.photonsDetectedb2c0 +  counterData.rollover.photonsDetectedb2c1 - counterData.rollover.photons_b2;
    
    %photons lost in processing
    counterData.rollover.photons_lost_b1 = diff(counterData.rollover.raw.beam1_photons_lost([1:end end]));
    counterData.rollover.photons_lost_b2 = diff(counterData.rollover.raw.beam2_photons_lost([1:end end]));
    
    %tags lost in processing
    counterData.rollover.tags_lost_b1 = diff(counterData.rollover.raw.beam1_tags_lost([1:end end]));
    counterData.rollover.tags_lost_b2 = diff(counterData.rollover.raw.beam2_tags_lost([1:end end]));

    %check for issues that can be diagnosed just from this file alone
    err.photonsLost = any(raw.beam1_photons_lost | raw.beam2_photons_lost);
    err.tagsLost = any(raw.beam1_tags_lost | raw.beam2_tags_lost);
    err.dataNotReceived = max(raw.dma_loop_bytes_sent_to_dma_fifo) > max(raw.elements_received_from_fpga)*2;
    err.countersSkipped = any(diff(counterData.rollover.cnum) ~= 1);
    counterData.errors = err;
    counterData.hasErrors = (err.photonsLost || err.tagsLost || err.dataNotReceived);
    
    fn = fieldnames(counterData);
    
    
    explanations = {'ms_timer - labview timer when this update was made'                                          ;
     'elements_received_from_fpga - number of u16s read to the time of this update'                     ;
     'pmt_loop_counter_-_160_MHz - counter value, this value and all following update every 2^28 ticks'                       ;
     'beam1_photon_u16s_in_to_pmt_volumetric - photons input to volumetric loop'           ;
     'beam1(/2)_tag_triggers_into_pmt_volumetric - tag triggers input to volumetric loop'            ;
     'beam1(/2)_total_tag_and_photons_sent_to_volumetric - total tags and photons that come out of interleave in encode photons'    ;
     'beam1(/2)_tags_lost - tags lost into fifo in encode photons'                                   ;
     'beam1(/2)_photons_lost - photons lost into fifo in encode photons'                                ;
     'beam1_u16s_to_interleave,beam2_u16s_to_interleave,total_beam1beam2_u16s_interleaved_into_photon_fifo - from interleave in encode photons dual beam'                          ;
     'total_u16s_sent_to_encode_timers_and_aod_signals - number of u16s out of counter & loc coding in volumetric loop dual channel'  ;
     'total_u16s_sent_to_encode_photons_and_tag_triggers - number of u16s out of photons FIFO in volumetric loop dual channel';
     'bytes_into_send_volumetric_as_u16 - bytes into send volumetric in volumetric loop dual channel'                 ;
     'bytes_into_transmission_fifo_as_u64 - bytes into FIFO in send volumetric via fifo (in volumetric loop dual channel)'               ;
     'b1(/2)ch0(/1)_num_pmt_signals_integrated - counts signals received from PMT independent of volumetric loop'                  ;
     'dma_loop_counter_-_100_MHz - counter inside 100 MHz dram access loop'                        ;
     'dma_loop_bytes_received_from_pmt_loop - bytes received by dram buffer loop'             ;
     'dma_loop_bytes_sent_to_dma_fifo - bytes out of dram buffer loop = bytes sent to computer'                   ;
     'consistency - ultimately dma_loop_bytes should equal elements_received*2, but because asynchronous will be ahead or behind until streaming stops' ;
      'consistency - elements_received*2 should equal file size' ;
      'consistency - difference between counter rollovers of bytes_into_send_volumetric_as_u16 should be same as difference in file location of coutner rollovers' ;
      'consistency - difference between #tag triggers, #photons between rollovers should be the same on fgpa and in file'; 
     
     };
end

