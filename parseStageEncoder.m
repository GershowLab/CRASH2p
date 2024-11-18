function enc = parseStageEncoder(encfname, eventfname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%'request_number	time_stamp	x_enc	y_enc	z_enc	f_enc	raw_buffer_hex'
j = 0;
fid = fopen(encfname);
fgetl(fid); %discard header line
try
    while(true)
        str = fgetl(fid);
        if (~ischar(str) || isempty(str))
            break;
        end
        j = j+1;
        [rn,~,~,ni] = sscanf(str, '%ld',1);
        enc.requestnum(j) = double(rn);
        str = str(ni:end);
        [mst,~,~,ni] = sscanf(str, '%ld',1);
        enc.mstime(j) = double(mst);
        str = str(ni:end);
        [elc,~,~,ni] = sscanf(str, '%ld', [4 1]);
        enc.eloc(:,j) = double(elc);
        enc.str{j} = str(ni:end);
        
        enc.hex{j} = sscanf(enc.str{j},'%2x');
    end
catch me
    fclose(fid);
    rethrow(me);
end

if (j == 0)
    warning ('no encoders in stage file');
    enc = [];
    return;
end

if (nargin < 2)
    eventfname = strrep(encfname, '_stage_enc.txt', '_event_log.txt');
end

ll = cellfun(@length, enc.hex);
enc.parsevalid = (~any(enc.eloc >= 2^31 -1, 1)) & ll == 21 & all(abs(enc.eloc) ~= 165); %fix weird error in stage code

for qq = 1:10
    validinds = find(enc.parsevalid);
    eloc = enc.eloc(:,validinds);
    valid = true(size(validinds));
    %discard points where encoder read is more than 1000 counts from the median
    %of its neighborhood 
    deltaencmax = 1000; 
    invloc = [];
    for j = 1:3
        valid = valid & (abs(eloc(j,:) - medfilt1(eloc(j,:),5)) < deltaencmax);
        if (any(~valid))
            invloc = [invloc mode(eloc(j,~valid))]; %#ok<AGROW>
        end
    end
    for j = 1:length(invloc)
        valid = valid & all(eloc ~= invloc(j));
    end
    enc.parsevalid(validinds) = valid;
    if (all(valid))
        break;
    end
end
    

STAGE_CODE = 2;
data = importdata(eventfname);
%{'ms_timer'  'code'  'data'  'timer'}
inds = data.data(:,2) == STAGE_CODE;
mstimer = data.data(inds,1);
reqnum = data.data(inds,3);
fpgatimer = data.data(inds, 4);

TIME_RANGE = 2000; %ms -- assume request was sent, received, processed, and reported back within 2 seconds;

enc.hasmatch = false(size(enc.requestnum));
enc.fpgatimer = NaN(size(enc.requestnum));
for j = 1:length(enc.requestnum)
    ind = find(enc.requestnum(j) == reqnum & abs(enc.mstime(j) - mstimer) < TIME_RANGE);
    if (~isempty(ind))
        [~,I] = min(abs(enc.mstime(j) - mstimer(ind)));
        enc.fpgatimer(j) = fpgatimer(ind(I));
        enc.hasmatch(j) = true;
    end
end

tcent = mean(enc.mstime(enc.hasmatch & enc.parsevalid));
tscale = max(enc.mstime(enc.hasmatch & enc.parsevalid)) - min(enc.mstime(enc.hasmatch & enc.parsevalid));

fcent = mean(enc.fpgatimer(enc.hasmatch & enc.parsevalid));
fscale = max(enc.fpgatimer(enc.hasmatch & enc.parsevalid)) - min(enc.fpgatimer(enc.hasmatch & enc.parsevalid));

p = polyfit((enc.mstime(enc.hasmatch & enc.parsevalid) - tcent)/tscale, (enc.fpgatimer(enc.hasmatch & enc.parsevalid)-fcent)/fscale, 1);
enc.fpgatimer(~enc.hasmatch) = polyval(p, (enc.mstime(~enc.hasmatch)-tcent)/tscale)*fscale + fcent;

enc.loc = NaN(3, length(enc.mstime));
XENC = 0.022;
YENC = XENC;
ZENC = 0.0055;
scale = repmat([XENC;YENC;ZENC], [1 length(enc.mstime)]);
enc.loc(:,enc.parsevalid) = scale(:,enc.parsevalid).*enc.eloc(1:3,enc.parsevalid);
enc.loc(:,~enc.parsevalid) = interp1(enc.mstime(enc.parsevalid), enc.loc(:,enc.parsevalid)', enc.mstime(~enc.parsevalid), 'linear', NaN)'; 
 enc.valid = all(isfinite(enc.loc));
% 
% try
%     stagefname = strrep(encfname, '_stage_enc.txt', '_stage.txt');
%     data = importdata(stagefname);
%     ch = matlab.lang.makeValidName(data.colheaders);
%     for j = 1:length(ch)
%         enc.stagefile.(ch{j}) = data.data(:,j);
%     end
%     
%     mstime = enc.stagefile.ms_timer;
%     %fpgtatime = polyval(p, (mstime-tcent)/tscale)*fscale + fcent;
%     stageloc = [enc.stagefile.x enc.stagefile.y enc.stagefile.z];
% 
%     enc.loc(:,~enc.parsevalid) = interp1(mstime, stageloc, enc.mstime(~enc.parsevalid), 'linear', NaN)'; 
%     enc.valid = all(isfinite(enc.loc));
% catch me
%     warning ('backup readings from stage.txt failed');
%     disp(me.getReport());
%       enc.valid = enc.parsevalid;% & enc.hasmatch;
% end
   
