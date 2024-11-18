classdef StreamingDataParseError < uint32
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    enumeration
        none (0), invalid_checksums (1), aod_in_timer (2), timer_in_aod (4), unexpected_rollover (8), missing_timers (16), aods_do_not_match (32), tags_do_not_match (64), 
    end
end

