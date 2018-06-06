function [mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
    movtm] = CellsortPCA(experiment, flims, nPCs, dsamp, outputdir, badframes, pw, ph, movieType)
% [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs, dsamp, outputdir, badframes)
%
% CELLSORT
% Read TIFF movie data and perform singular-value decomposition (SVD)
% dimensional reduction.
%
% Inputs:
%   experiment - experiment structure
%   flims - 2-element vector specifying the endpoints of the range of
%   frames to be analyzed. If empty, default is to analyze all movie
%   frames.
%   nPCs - number of principal components to be returned
%   dsamp - optional downsampling factor. If scalar, specifies temporal
%   downsampling factor. If two-element vector, entries specify temporal
%   and spatial downsampling, respectively.
%   outputdir - directory in which to store output .mat files (relative to the experiment folder)
%   badframes - optional list of indices of movie frames to be excluded
%   from analysis
%   pw - 2-element vector with first and last pixel width coordiantes to use (leave empty to use the whole image)
%   ph - 2-element vector with first and last pixel height coordiantes to use (leave empty to use the whole image)
%   movieType - original/glia
%
% Outputs:
%   mixedsig - N x T matrix of N temporal signal mixtures sampled at T
%   points.
%   mixedfilters - N x X x Y array of N spatial signal mixtures sampled at
%   X x Y spatial points.
%   CovEvals - largest eigenvalues of the covariance matrix
%   covtrace - trace of covariance matrix, corresponding to the sum of all
%   eigenvalues (not just the largest few)
%   movm - average of all movie time frames at each pixel
%   movtm - average of all movie pixels at each time frame, after
%   normalizing each pixel deltaF/F
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

%
% June-August 2016 (Cellsort 1.41): Fixed an inconsistency between create_tcov and
% creat_xcov. In the current version, both methods agree with the previous
% version of create_tcov.
%

% 2017 - NETCAL Adaptation - Javier Orlandi


tic
fprintf('-------------- CellsortPCA %s: %s -------------- \n', date, experiment.name)

%-----------------------
% Check inputs
% if isempty(dir(fn))
%     error('Invalid input file name.')
% end
if (nargin<2)||(isempty(flims))
    nt_full = experiment.numFrames;
    flims = [1,nt_full];
end
if nargin<6
    badframes = [];
end

useframes = setdiff((flims(1):flims(2)), badframes);
nt = length(useframes);

if nargin<3 || isempty(nPCs)
    nPCs = min(150, nt);
end
if nargin<4 || isempty(dsamp)
    dsamp = [1,1];
end
if nargin<5 || isempty(outputdir)
    outputdir = [experiment.folder, filesep 'cellsort_preprocessed_data', filesep];
end
if isempty(dir(outputdir))
    mkdir([experiment.folder, filesep 'cellsort_preprocessed_data', filesep]);
end
if outputdir(end)~= filesep
    outputdir = [outputdir, filesep];
end

%[~, fname] = fileparts(fn);
fname = experiment.name;
if isempty(badframes)
    fnmat = [outputdir, fname, '_',num2str(flims(1)),',',num2str(flims(2)), '_', date,'.mat'];
else
    fnmat = [outputdir, fname, '_',num2str(flims(1)),',',num2str(flims(2)),'_selframes_', date,'.mat'];
end
if ~isempty(dir(fnmat))
    fprintf('CELLSORT: Movie %s already processed;', ...
        experiment.name)
    forceload = input(' Re-load data? [0-no/1-yes] ');
    if isempty(forceload) || forceload==0
        load(fnmat)
        return
    end
end

if(nargin < 7)
  pw = 1:experiment.width;
else
  pw = pw(1):pw(2);
end
if(nargin < 8)
  ph = 1:experiment.height;
else
  ph = ph(1):ph(2);
end
pixw = numel(pw);
pixh = numel(ph);
npix = pixw*pixh;

fprintf('   %d pixels x %d time frames;', npix, nt)

% Create covariance matrix
if nt < npix 
    fprintf(' using temporal covariance matrix.\n')
    [covmat, mov, movm, movtm] = create_tcov(experiment, pw, ph, useframes, nt, dsamp, movieType);
else
    fprintf(' using spatial covariance matrix.\n')
    [covmat, mov, movm, movtm] = create_xcov(experiment, pw, ph, useframes, nt, dsamp, movieType);
end

covtrace = trace(covmat) / npix;
movm = reshape(movm, pixw, pixh);

if nt < npix 
    % Perform SVD on temporal covariance
    [mixedsig, CovEvals] = cellsort_svd(covmat, nPCs, nt, npix);

    % Load the other set of principal components
    [mixedfilters] = reload_moviedata(pixw*pixh, mov, mixedsig, CovEvals);
else
    % Perform SVD on spatial components
    [mixedfilters, CovEvals] = cellsort_svd(covmat, nPCs, nt, npix);
    mixedfilters = mixedfilters' * npix;
    
    % Load the other set of principal components
    [mixedsig] = reload_moviedata(nt, mov', mixedfilters', CovEvals);
    mixedsig = mixedsig' / npix^2;
end
%mixedfilters = reshape(mixedfilters, pixw,pixh,nPCs);
mixedfilters = reshape(mixedfilters, pixh,pixw,nPCs);

%------------
% Save the output data
save(fnmat,'mixedfilters','CovEvals','mixedsig', ...
    'movm','movtm','covtrace')
fprintf(' CellsortPCA: saving data and exiting; ')
toc

    function [covmat, mov, movm, movtm] = create_xcov(experiment, pw, ph, useframes, nt, dsamp, movieType)
        %-----------------------
        % Load movie data to compute the spatial covariance matrix
        pixw1 = numel(pw);
        pixh1 = numel(ph);
        npix1 = pixw1*pixh1;
        % Check if we are not using the whole frame
        if(pixw1 ~= experiment.width || pixh1 ~= experiment.height)
          % Inverted because of HIS structure
          [rf,cf]= find(ones(experiment.height, experiment.width));
          pixelList = find(cf >= min(ph) & cf <= max(ph) & rf >= min(pw) & rf <= max(pw));
        else
          pixelList = [];
        end
        % Load the video Stream
        switch movieType
          case 'glia'
            fID = fopen([experiment.folder 'data' filesep experiment.name '_gliaAverageMovieDataFile.dat'], 'r');
            fseek(fID, 0, 'bof'); % Just in case
            width = size(experiment.gliaAverageFrame, 2);
            height = size(experiment.gliaAverageFrame, 2);
          otherwise
            [fID, experiment] = openVideoStream(experiment);
        end
        
        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end

        if (dsamp_space==1)
            mov = zeros(pixh1, pixw1, nt);
            for jjind=1:length(useframes)
                jj = useframes(jjind);
                %mov(:,:,jjind) = imread(fn,jj);
                switch movieType
                  case 'glia'
                    fseek(fID, prod(width*height)*8*(jj-1), 'bof'); % the 8 refers to BYTES
                    tmpData = fread(fID, [size(experiment.gliaAverageFrame, 2) size(experiment.gliaAverageFrame, 1)], 'double'); % Sequential read
                    curFrame = tmpData(pixelList);
                  otherwise
                    curFrame = getFrame(experiment, jj, fID, pixelList);
                end
                if(~isempty(pixelList))
                  %curFrame = reshape(curFrame, [pixh1, pixw1]);
                  % Inverted because of HIS ordering
                  curFrame = reshape(curFrame, [pixw1, pixh1])';
                end
                mov(:, :, jjind) = curFrame;
                if mod(jjind,500)==1
%                   figure;
%                   imagesc(curFrame);
%                   colormap gray;
%                   colorbar;
%                   [m, M] = autoLevelsFIJI(curFrame, 16, true);
%                   caxis([m, M]);
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        else
          switch movieType
            case 'glia'
              fseek(fID, prod(width*height)*8*(0), 'bof'); % the 8 refers to BYTES
              tmpData = fread(fID, [size(experiment.gliaAverageFrame, 2) size(experiment.gliaAverageFrame, 1)], 'double'); % Sequential read
              curFrame = tmpData(pixelList);
            otherwise
              curFrame = getFrame(experiment, 1, fID, pixelList);
          end
          %curFrame = getFrame(experiment, 1, fID, pixelList);
          if(~isempty(pixelList))
            curFrame = reshape(curFrame, [pixh1, pixw1]);
          end
          [pixw_dsamp,pixh_dsamp] = size(imresize(curFrame, 1/dsamp_space, 'bilinear' ));
          mov = zeros(pixw_dsamp, pixh_dsamp, nt);
          for jjind=1:length(useframes)
              jj = useframes(jjind);
              %mov(:,:,jjind) = imresize( imread(fn,jj), 1/dsamp_space, 'bilinear' );
              switch movieType
                case 'glia'
                  fseek(fID, prod(width*height)*8*(jj-1), 'bof'); % the 8 refers to BYTES
                  tmpData = fread(fID, [size(experiment.gliaAverageFrame, 2) size(experiment.gliaAverageFrame, 1)], 'double'); % Sequential read
                  curFrame = tmpData(pixelList);
                otherwise
                  curFrame = getFrame(experiment, jj, fID, pixelList);
              end
              if(~isempty(pixelList))
                curFrame = reshape(curFrame, [pixh1, pixw1]);
              end
              mov(:, :, jjind) = imresize(curFrame, 1/dsamp_space, 'bilinear');
              if mod(jjind,500)==1
                  fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                  toc
              end
          end
        end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix1, nt);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0);
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1; % Compute Delta F/F
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp,1)/dsamp, 1, mov, [], 2);
            mov = downsample(mov', dsamp)';
        end

        movtm = mean(mov,1); % Average over space
        mov = mov - ones(size(mov,1),1)*movtm; 

        covmat = (mov*mov')/size(mov,2);
        covmat = covmat * size(mov,2)/size(mov,1); % Rescale to gree with create_tcov
        toc
        closeVideoStream(fID);
    end

    function [covmat, mov, movm, movtm] = create_tcov(experiment, pw, ph, useframes, nt, dsamp, movieType)
        %-----------------------
        % Load movie data to compute the temporal covariance matrix
        pixw1 = numel(pw);
        pixh1 = numel(ph);
        npix1 = pixw1*pixh1;
        % Check if we are not using the whole frame
        if(pixw1 ~= experiment.width || pixh1 ~= experiment.height)
          [c, r] = meshgrid(pw, ph);
          pixelList = sub2ind([experiment.height, experiment.width], r(:), c(:));
        else
          pixelList = [];
        end
        % Load the video Stream
        switch movieType
          case 'glia'
            width = size(experiment.gliaAverageFrame, 2);
            height = size(experiment.gliaAverageFrame, 2);
            fID = fopen([experiment.folder 'data' filesep experiment.name '_gliaAverageMovieDataFile.dat'], 'r');
            fseek(fID, 0, 'bof'); % Just in case
          otherwise
            [fID, experiment] = openVideoStream(experiment);
        end
        
        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end

        if (dsamp_space==1)
            %mov = zeros(pixw1, pixh1, nt);
            mov = zeros(pixh1, pixw1, nt);
            for jjind=1:length(useframes)
                jj = useframes(jjind);
                %mov(:,:,jjind) = imread(fn,jj);
                switch movieType
                  case 'glia'
                    fseek(fID, prod(width*height)*8*(jj-1), 'bof'); % the 8 refers to BYTES
                    tmpData = fread(fID, [size(experiment.gliaAverageFrame, 2) size(experiment.gliaAverageFrame, 1)], 'double'); % Sequential read
                    curFrame = tmpData(pixelList);
                  otherwise
                    curFrame = getFrame(experiment, jj, fID, pixelList);
                end
                if(~isempty(pixelList))
                  curFrame = reshape(curFrame, [pixh1, pixw1]);
                end
                mov(:, :, jjind) = curFrame;
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        else
          switch movieType
            case 'glia'
              fseek(fID, prod(width*height)*8*(0), 'bof'); % the 8 refers to BYTES
              tmpData = fread(fID, [size(experiment.gliaAverageFrame, 2) size(experiment.gliaAverageFrame, 1)], 'double'); % Sequential read
              curFrame = tmpData(pixelList);
            otherwise
              curFrame = getFrame(experiment, 1, fID, pixelList);
          end
          if(~isempty(pixelList))
            curFrame = reshape(curFrame, [pixh1, pixw1]);
          end
          [pixw_dsamp,pixh_dsamp] = size(imresize(curFrame, 1/dsamp_space, 'bilinear' ));
          mov = zeros(pixw_dsamp, pixh_dsamp, nt);
          for jjind=1:length(useframes)
              jj = useframes(jjind);
              %mov(:,:,jjind) = imresize( imread(fn,jj), 1/dsamp_space, 'bilinear' );
              switch movieType
                case 'glia'
                  fseek(fID, prod(width*height)*8*(jj-1), 'bof'); % the 8 refers to BYTES
                  tmpData = fread(fID, [size(experiment.gliaAverageFrame, 2) size(experiment.gliaAverageFrame, 1)], 'double'); % Sequential read
                  curFrame = tmpData(pixelList);
                otherwise
                  curFrame = getFrame(experiment, jj, fID, pixelList);
              end
              if(~isempty(pixelList))
                curFrame = reshape(curFrame, [pixh1, pixw1]);
              end
              mov(:, :, jjind) = imresize(curFrame, 1/dsamp_space, 'bilinear');
              if mod(jjind,500)==1
                  fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                  toc
              end
          end
        end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix1, nt);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0); % Avoid dividing by zero
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1;
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp,1)/dsamp, 1, mov, [], 2);
            mov = downsample(mov', dsamp)';
        end

        c1 = (mov'*mov)/npix1;
        movtm = mean(mov,1); % Average over space
        covmat = c1 - movtm'*movtm;
        clear c1
        closeVideoStream(fID);
    end

    function [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix1)
        %-----------------------
        % Perform SVD

        covtrace1 = trace(covmat) / npix1;

        opts.disp = 0;
        opts.issym = 'true';
        if nPCs<size(covmat,1)
            [mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM', opts);  % pca_mixedsig are the temporal signals, mixedsig
        else
            [mixedsig, CovEvals] = eig(covmat);
            CovEvals = diag( sort(diag(CovEvals), 'descend'));
            nPCs = size(CovEvals,1);
        end
        CovEvals = diag(CovEvals);
        if nnz(CovEvals<=0)
            nPCs = nPCs - nnz(CovEvals<=0);
            fprintf(['Throwing out ',num2str(nnz(CovEvals<0)),' negative eigenvalues; new # of PCs = ',num2str(nPCs),'. \n']);
            mixedsig = mixedsig(:,CovEvals>0);
            CovEvals = CovEvals(CovEvals>0);
        end

        mixedsig = mixedsig' * nt;
        CovEvals = CovEvals / npix1;

        percentvar = 100*sum(CovEvals)/covtrace1;
        fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])
    end

    function [mixedfilters] = reload_moviedata(npix1, mov, mixedsig, CovEvals)
        %-----------------------
        % Re-load movie data
        nPCs1 = size(mixedsig,1);

        Sinv = inv(diag(CovEvals.^(1/2)));

        movtm1 = mean(mov,1); % Average over space
        movuse = mov - ones(npix1,1) * movtm1;
        mixedfilters = reshape(movuse * mixedsig' * Sinv, npix1, nPCs1);
        %mixedfilters = reshape(movuse * mixedsig * Sinv, npix1, nPCs1);
    end
end