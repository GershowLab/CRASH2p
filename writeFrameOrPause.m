function writeFrameOrPause(vidObj, avirect, fps, nframes, flushBuffer)
%function writeFrameOrPause(vidObj, avirect, fps, nframes, flushBuffer)
%batch processes frames to speed up writing
%flushBuffer writes current buffer to disk, without adding any more frames
%to the buffer
%this should be redone as an object bufferedVideoWriter
%
%https://www.mathworks.com/matlabcentral/answers/1929240-how-to-speed-up-a-script-writing-frames-using-a-videowriter-object
%
%call with no arguments to reset state

existsAndDefault('fps', 25);
existsAndDefault('nframes', 1);
existsAndDefault('avirect', []);
existsAndDefault('flushBuffer', false);

stackDepth = 200;
persistent frameStack;
persistent ind;

if (nargin <= 0)
    frameStack = [];
    ind = 0;
    return;
end

if (flushBuffer && ~isempty(vidObj) && ~isempty(frameStack))
    writeVideo(vidObj,frameStack(1:ind));
    return;
end


for j = 1:nframes
    if (~isempty(vidObj))
        
        if (isempty(frameStack))
            if (isempty(avirect))
                currFrame = getframe(gcf);
            else
                currFrame = getframe(gcf,avirect);
            end
            frameStack = repmat(currFrame, stackDepth, 1);
            ind = 1;
        else
            if (isempty(avirect))
                frameStack(ind) =getframe(gcf);
            else
                frameStack(ind) = getframe(gcf,avirect);
            end
        end
        ind = ind + 1;
       if (ind > stackDepth)
           writeVideo(vidObj,frameStack);
           ind = 1;
       end
       
    else
        pause(1/fps);
    end
end