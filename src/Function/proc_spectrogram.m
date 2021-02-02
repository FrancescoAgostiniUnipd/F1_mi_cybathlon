function [features, f] = proc_spectrogram(data, wlenght, wshift, pshift, samplerate, mlength)
% [features, f] = proc_spectrogram(data, wlenght, wshift, pshift, samplerate [, mlength])
%
% The function computes the spectrogram on the real data.
%
% Input arguments:
%   - data              Data matrix [samples x channels]
%   - wlength           Window's lenght to be used to segment data and
%                       compute the spectrogram                             [in seconds]
%   - wshift            Shift of the external window (e.g., frame size)     [in seconds]
%   - pshift            Shift of the internal psd windows                   [in seconds]
%   - samplerate        Samplerate of the data
%   - [mlength]         Optional length of the external windows to compute
%                       the moving average.                                 [in seconds] 
%                       By default the length of the moving average window
%                       is set to 1 second. To not compute the moving
%                       average, empty argument can be provided.
%
% Output arguments:
%   - features          Output of the spectrogram in the format: 
%                       [windows x frequencies x channels]. Number of
%                       windows (segments) is computed according to the
%                       following formula: 
%                       nsegments = fix((NX-NOVERLAP)/(length(WINDOW)-NOVERLAP))
%                       where NX is the total number of samples, NOVERLAP
%                       the number of overlapping samples for each segment
%                       and length(WINDOW) the number of samples in each
%                       segment. 
%                       Number of frequencies is computed according to the
%                       NFFT. nfrequencies is equal to (NFFT/2+1) if NFFT 
%                       is even, and (NFFT+1)/2 if NFFT is odd. NFFT is the
%                       maximum between 256 and the next power of 2 greater
%                       than the length(WINDOW).
%   - f                 Vectore with the computed frequencies
%
% SEE ALSO: spectrogram, nextpow2
    
    % Data informations
    nsamples  = size(data, 1);
    nchannels = size(data, 2);

    % Useful params for PSD extraction with the fast algorithm
    psdshift = pshift*samplerate;
    winshift = wshift*samplerate;

    if((mod(psdshift,winshift) ~=0) && (mod(winshift,psdshift) ~=0))
        warning('chk:par', '[proc_spectrogram] The fast PSD method cannot be applied with the current settings!');
        error('chk:par', '[proc_spectrogram] The internal welch window shift must be a multiple of the overall feature window shift (or vice versa)!');
    end

    % Create arguments for spectrogram
    spec_win = wlenght*samplerate;
    
    % Careful here: The overlapping depends on whether the winshift or the
    % psdshift is smaller. Some calculated internal windows will be redundant,
    % but the speed is much faster anyway

    if(psdshift <= winshift)
        spec_ovl = spec_win - psdshift;
    else
        spec_ovl = spec_win - winshift;
    end

    % Calculate all the internal PSD windows
    nsegments = fix((nsamples-spec_ovl)/(spec_win-spec_ovl));          % From spectrogram's help page
    nfft      = max(256, nextpow2(spec_win));
    if(mod(nfft, 2) == 0)
        nfreqs = (nfft/2) + 1;
    else
        nfreqs = (nfft+1)/2;
    end
    
    psd = zeros(nfreqs, nsegments, nchannels);
    for chId = 1:nchannels
        [~,f,~,psd(:,:,chId)] = spectrogram(data(:,chId), spec_win, spec_ovl, [], samplerate);
    end
    
    if nargin == 5
        mlength = 1;
    end
    
    if isempty(mlength) == false
        % Setup moving average filter parameters
        mavg_a = 1;
        if(winshift >= psdshift)
            % Case where internal windows are shifted according to psdshift
            mavgsize  = ((mlength*samplerate)/psdshift) - 1;   
            mavg_b    = (1/mavgsize)*ones(1,mavgsize);
            mavg_step = winshift/psdshift;
        else
            % Case where internal windows are shifted according to winshift
            mavgsize  = ((mlength*samplerate)/winshift) - (psdshift/winshift);   
            mavg_b    = zeros(1,mavgsize);
            mavg_b(1:psdshift/winshift:end-1) = 1;
            mavg_b    = mavg_b/sum(mavg_b);
            mavg_step = 1;
        end
        
        startindex = find(mavg_b~=0, 1, 'last');

        features = filter(mavg_b,mavg_a,psd,[],2);
        features = permute(features, [2 1 3]);

        % Get rid of initial filter byproducts
        features = features(startindex:end,:,:);

        % In case of psdshift, there will be redundant windows. Remove them
        if(mavg_step > 1)
           features = features(1:mavg_step:end,:,:);
        end
    else
        features = psd;
        features = permute(features, [2 1 3]);
    end
end