function wPOS = proc_pos2win(POS, wshift, direction, wlength)
% wPOS = proc_pos2win(POS, wshift, edge, wlength)
%
% The function converts sample position (POS) to window (psd) position for 
% a given  window shift (wshift). direction argument refers to the way the
% conversion is done (according to the past or future (psd) window. If 
% backward direction is required, then also the length of the window must be
% provided.
%
% Input:
%   - POS           Vector with sample positions
%   - wshift        Window shift                                            [samples]
%   - direction     Either 'forward' or 'backward'. Forward means that the
%                   conversion assumes to have data in the future window
%                   (with respect to the position). Backward means that the
%                   position refers to the past data.
%   - wlength       Required argument if rising edge is requested. It
%                   corresponds to the number of samples of the window.     [samples]
%
% Output:
%   - wPOS          Window position
% 
% SEE ALSO: proc_spectrogram

    backward = false;
    if strcmpi(direction, 'forward')
        wlength = [];
    elseif strcmpi(direction, 'backward')
        backward = true;
        if nargin == 3
            error('chk:arg', 'backward direction option requires to provide wlength');
        end
    else
        error('chk:arg', 'Direction not recognized: only forward and backward are allowed');
    end
    
    wPOS = floor(POS/wshift) + 1;
    
    if backward == true
        wPOS = wPOS - (floor(wlength/wshift));
    end

end