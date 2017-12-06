function cell_sig = CellsortApplyFilter(experiment, ica_segments, flims, movm, subtractmean, pw, ph)
% cell_sig = CellsortApplyFilter(fn, ica_segments, flims, movm, subtractmean)
%
%CellsortApplyFilter
% Read in movie data and output signals corresponding to specified spatial
% filters
%
% Inputs:
%   experiment - experiment structure
%     ica_segments - nIC x X matrix of ICA spatial filters
%     flims - optional two-element vector of frame limits to be read
%     movm - mean fluorescence image
%     subtractmean - boolean specifying whether or not to subtract the mean
%     fluorescence of each time frame
%   pw - 2-element vector with first and last pixel width coordiantes to use (leave empty to use the whole image)
%   ph - 2-element vector with first and last pixel height coordiantes to use (leave empty to use the whole image)
% Outputs:
%     cell_sig - cellular signals
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%
%
% 2017 - NETCAL Adaptation - Javier Orlandi


if (nargin<3)||isempty(flims)
    nt = experiment.numFrames;
    flims = [1,nt];
else
    nt = diff(flims)+1;
end
if nargin<5
    subtractmean = 1;
end
if(nargin < 6)
  pw = 1:experiment.width;
else
  pw = pw(1):pw(2);
end
if(nargin < 7)
  ph = 1:experiment.height;
else
  ph = ph(1):ph(2);
end
pixw = numel(pw);
pixh = numel(ph);


if (nargin<4)||isempty(movm)
    movm = ones(pixh,pixw);
else
    movm = double(movm);
end
movm(movm==0) = 1; % Just in case there are black areas in the average image
k=0;

cell_sig = zeros(size(ica_segments,1), nt);
ica_segments = reshape(ica_segments, [], pixw*pixh);

if(pixw ~= experiment.width || pixh ~= experiment.height)
  % Inverted because of HIS structure
  [rf,cf]= find(ones(experiment.height, experiment.width));
  pixelList = find(cf >= min(ph) & cf <= max(ph) & rf >= min(pw) & rf <= max(pw));
else
  pixelList = [];
end
% Load the video Stream
[fID, experiment] = openVideoStream(experiment);
        
fprintf('Loading %5g frames for %d ROIs.\n', nt, size(ica_segments,1))
tic
while k<nt
    ntcurr = min(500, nt-k);
    mov = zeros(pixh, pixw, ntcurr);
    for j=1:ntcurr
        %movcurr = imread(fn, j+k+flims(1)-1);
        curFrame = getFrame(experiment, j+k+flims(1)-1, fID, pixelList);
        if(~isempty(pixelList))
          % Inverted because of HIS ordering
          curFrame = reshape(curFrame, [pixw, pixh])';
        end

        mov(:,:,j) = curFrame;
    end
    mov = (mov ./ repmat(movm, [1,1,ntcurr])) - 1; % Normalize by background and subtract mean
    
    if subtractmean
        % Subtract the mean of each frame
        movtm = mean(mean(mov,1),2);
        mov = mov - repmat(movtm,[pixh,pixw,1]);
    end
    
    mov = reshape(mov, pixw*pixh, ntcurr);
    cell_sig(:, k+[1:ntcurr]) = ica_segments*mov;
    
    k=k+ntcurr;
    fprintf('Loaded %3.0f frames; ', k)
    toc
end
closeVideoStream(fID);
