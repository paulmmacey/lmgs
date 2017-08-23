function cspm_lmgs(Pin,overwrite,prefix,display,GM)
    % CSPM_LMGS - Detrend fMRI image series using LMGS method.
    % CSPM_LMGS(P, OVERWRITE, PREFIX, DISPLAY, GRANDMEAN )
    %       P is array of images returned by spm, e.g. (SPM2):
    %           P = spm_get(Inf, '.img',{'Please select images for detrending'});
    %           You will be prompted for the files if the argument is omitted.
    %           To pass multiple sessions, set P as a cell array with each cell the
    %           files for one session, i.e. (SPM5):
    %           P{1} = = spm_select([3 Inf],'image','Select images for detrending (1st session)');
    %           P{2} = = spm_select([3 Inf],'image','Select images for detrending (2nd session)');
    %           P{3} = = spm_select([3 Inf],'image','Select images for detrending (3rd session)');
    %       OVERWRITE (optional) - defaults to true; if false and detrended
    %           images exist, detrending is skipped.
    %       PREFIX (optional) - prefix for detrended files; defaults to "d".
    %       DISPLAY (optional) - if true, global trends of raw and detrended
    %           data are plotted; defaults to false.
    %       GRANDMEAN (optional) - scale series to this mean (e.g., 100); default 0 (no scaling);
    %
    % The procedure applies the detrending to the time-series of files P, and writes
    % detrended files prepended with PREFIX (defaults to "d").
    %
    % To run in interactive mode, type
    %           >>cspm_lmgs
    % Interactive mode defaults to overwrite TRUE and prefix "d".
    %
    % If at least one input argument is provided (P), the routine will run in batch
    % mode and use the defaults of any arguments not passed.
    %
    %
    % LMGS (Linear Model of Global Signal) detrending is described in:
    % P.M. Macey, K.E. Macey, R. Kumar, and R.M. Harper. A method for removal of
    % global effects from fMRI time series. NeuroImage 22 (2004) 360-366.
    %
    % Works with SPM5+ and SPM2.
    %

    % cspm_lmgs.m	1.4 Paul Macey 2017-08-23
    % Thanks to Marko Wilke for the batch mode and grand mean scaling suggestions.

    disp('========================================================================')
    disp('LMGS-Detrending')
    disp('  ')

    % Check inputs
    % Images
    if nargin == 0
        overwrite = 1;
        prefix = 'd';
        num_sessions = spm_input('Number of sessions to detrend', 1, 'n', 1);
        Pall = cell(1,num_sessions);
        for i = 1:num_sessions
            if strcmp(spm('ver'),'SPM99') || strcmp(spm('ver'),'SPM2')
                Pall{i} = spm_get([3 Inf], '.img',['Select images for detrending (session ',...
                    int2str(i),')']);
            else
                Pall{i} = spm_select([3 Inf],'image',['Select images for detrending (session ',...
                    int2str(i),')']);
            end

            % Couple of checks
            if isempty(Pall{i}), return, end
            if length(Pall{i}) < 3
                msgbox('Detrending requires a time-series of images (multiple volumes)',...
                    'LMGS',...
                    'warn')
                i = i-1; %#ok<FXSET>
                if num_sessions == 1
                    return
                end
            end

        end
        display = spm_input('Show graphical results ?','+1','b',(' Yes| No'),[1 0],0);
        disp('Grand Mean scaling sets session mean to GM')
        GM = spm_input('Grand mean session scaling (0 = none)', '+1', 'e', 0);

    else
        % batch or silent mode
        if iscell(Pin)
            % Check case of cell array of filenames (i.e., one session)
            if min(size(Pin{1})) == 1
                Pall{1} = char(Pin);
            else
                Pall = Pin;
            end
        else
            % One session only
            Pall{1} = Pin;
        end
        num_sessions = length(Pall);
        % Flags, etc.
        if nargin < 2 || isempty(overwrite)
            overwrite = 1;
        end
        if nargin < 3 || isempty(prefix)
            prefix = 'd';
        else
            if ~ischar(prefix)
                msgbox('Invalid "Prefix" argument - need character',...
                    'LMGS','error')
                help cspm_lmgs
                return
            end
        end
        if nargin < 4
            display = 0;
        end
        if nargin < 5
            GM = 0;
        end
    end


    for i = 1:num_sessions
        if num_sessions > 1
            disp(['Detrending session ',int2str(i)])
        end
        detrend_onesession(Pall{i},overwrite,prefix,display,GM,i,num_sessions)
    end

    disp('LMGS-Detrending complete.')
    disp(' ')
end
% =========================================================================
function detrend_onesession(P,overwrite,prefix,display,GM,sn,num_sessions)

    % Number of scans
    nscan = size(P,1);


    % Check for 4D
    is4D = 0;
    [pth,nm,ext] = fileparts(P(1,:));
    I = find(ext==',',1,'first');
    if ~isempty(I)
        V = spm_vol([pth,filesep,nm,ext(1:I-1)]);
        if length(V) > 1
            is4D = 1;
            disp('4D file - assuming all images belong to same 4D file')
        end
    end

    % Finish now if overwrite flag is false and files exist
    if ~overwrite
        skip = 1;

        % Check for existing files - if 4D only have one file
        if is4D
            thiscount = 1;
        else
            thiscount = nscan;
        end
        for i = 1:thiscount

            [pth,nm,ext] = fileparts(P(i,:));

            % Detrended
            I = find(ext==',',1,'first');
            if ~isempty(I)
                ext = ext(1:I-1);
            end
            newfname = [pth,filesep,prefix,nm,ext];

            if ~exist(newfname,'file')
                skip = 0;
                break
            end
        end
        if skip, return, end
    end

    % Map volumes
    V = spm_vol(P);
    % Check similar dimensions, etc. (Taken from spm_imcalc_ui.)
    flagDimOK = true;
    if (strcmp(spm('ver'),'SPM99') || strcmp(spm('ver'),'SPM2')) 
        if any(any(diff(cat(1,V.dim),1,1),1)&[1,1,1,0]) 
            flagDimOK = false;
        end
    elseif any(any(diff(cat(1,V.dim),1,1),1))
        flagDimOK = false;        
    elseif any(any(any(diff(cat(3,V.mat),1,3),3)))
        flagDimOK = false;
    end
    if ~flagDimOK
        msgbox('Images don''t have same dimmensions, orientation and voxel size - detrending cancelled.',...
            'LMGS','error')
        return
    end

    % Number of voxels
    num_vox = prod(V(1).dim(1:3));

    % Initialize - create independent variable X for model
    % Global values
    gt = cspm_globaltrend(P,0,0,0)';
    % Remove offset
    X = detrend(gt,'constant');
    % Add Constant
    X = [X ones(nscan,1)];

    % Create output images
    % Preallocate
    Vout = repmat( V(1), 1, nscan );
    for i = 1:nscan
        Vnew = V(i);
        Vnew.dt(1) = spm_type('float32');
        [pth,nm,ext] = fileparts(Vnew.fname);
        I = find(ext==',',1,'first');
        if ~isempty(I)
            ext = ext(1:I-1);
        end
        if is4D
            str = ['3D',num2str(i,'%03d')];
            Vnew.n(1) = 1;
        else
            str = '';
        end
        Vnew.fname = [pth,filesep,prefix,nm,str,ext];
        if strcmp(spm('ver'),'SPM99')
            Vout(i) = spm_create_image(Vnew);
        else
            Vout(i) = spm_create_vol(Vnew);
        end
    end

    adjustcount = 0;
    num_vox_per_plane = num_vox / V(1).dim(3);
    dim1 = V(1).dim(1);
    dim2 = V(1).dim(2);
    dim3 = V(1).dim(3);
    dim12 = V(1).dim(1:2);
    % Loop through planes
    fprintf('\n%-60s','Detrending')
    for z = 1:dim3
        fprintf('%s%-60s',char(sprintf('\b')*ones(1,60)),...
            ['Detrending plane ',int2str(z),' / ',int2str(dim3)])

        % Get raw data for plane across all volumes.
        yn = zeros(dim1,dim2,nscan);
        for i = 1:nscan
            yn(:,:,i) = spm_slice_vol(V(i),spm_matrix([0 0 z]),dim12,0);
        end

        % Use 2-D array to speed processing.
        yadj = reshape(yn,num_vox_per_plane,nscan)';

        s = warning; warning off 
        % Perform detrending voxel-by-voxel in this plane
        for v = 1:num_vox_per_plane
            yvox = yadj(:,v);
            if any(yvox)
                b = regression(yvox,X);
                % The following would require the statistics toolbox:
                %                 [b,bint,r,rint,stats] = regress(yvox,X,0.05);
                % Test for positive effect - b(1);
                % Test for significance as p-value - - stats(3);
                % Test for correlation (R^2) - stats(1)                
                % Optional: skip voxels that do not meet thresholds
                %                 if b(1) <= 0 | stats(3) > 0.8 | stats(1) < 0.3
                %                     continue
                %                 end
                % Another option: remove only voxels that show positive correlations
                %                 if b(1) <= 0
                %                     continue
                %                 end
                
                yvoxmodel = b(1)*X(:,1);

                % Adjust time-series
                adjustcount = adjustcount + 1;
                yadj(:,v) = yvox - yvoxmodel;

            end
        end
        warning(s)
        yadj = reshape(yadj',dim1,dim2,nscan);

        % Write adjusted data
        for i = 1:nscan
            Vout(i) = spm_write_plane(Vout(i),yadj(:,:,i),z);
        end

    end     % End loop through planes
    fprintf('%s%-60s',char(sprintf('\b')*ones(1,60)),'Detrending')
    fprintf('%-60s\n','Done')

    disp(['Adjusted ', int2str(adjustcount),...
        ' non-zero voxels (',int2str(num_vox),' total).'])
    disp(' ')

    % Save file names if displaying or scaling (or if 4D, saving back to single 4D
    % file)
    if display || GM > 0 || is4D
        Pd = cell(nscan,1);
        for i = 1:nscan
            % Pin{i} = V(i).fname;
            Pd{i} = Vout(i).fname;
        end
    end

    % Close volumes, if SPM2
    if strcmp(spm('ver'),'SPM2')
        spm_close_vol(V);
        spm_close_vol(Vout);
    end

    % Grand mean scaling
    if GM > 0
        % Calculate required scale to set session mean to GM
        scale =  GM/(mean(cspm_globaltrend(Pd,0,0,0)));
        % Scale each volume by this amount
        for i = 1:nscan
            V = spm_vol(Pd{i});
            V.pinfo(1:2,:) = V.pinfo(1:2,:)*scale;
            V = spm_create_vol(V);
            if strcmp(spm('ver'),'SPM2')
                spm_close_vol(V);
            end
        end
    end

    % Show in percent change the pre- and post-detrending global timetrends
    if display

        % Get trends
        t1 = gt';   % Variable gt was calculated earlier.
        t2 = cspm_globaltrend(Pd,0,0,0)';

        b1 = mean(t1(1:end));
        b2 = mean(t2(1:end));

        % Percent change
        t1pc = 100*(t1 - b1)/b1;
        t2pc = 100*(t2 - b2)/b2;

        % For x-axis
        x = 1:length(gt);

        % Figure
        if num_sessions > 1
            titlestr = ['Effect of adjustment (session ',int2str(sn),')'];
        else
            titlestr = 'Effect of adjustment';
        end
        figure('NumberTitle','off','Name',titlestr)

        subplot(2,1,1)
        plot(x,t1,'b.-',x,t2,'r.-')
        title ('Raw values: Original and detrended images')
        legend('Raw (blue)','Adjusted (red)')

        subplot(2,1,2)
        plot(x,t1pc,'b.-',x,t2pc,'r.-')
        title ('Percent change: Original and detrended images (%)')
        legend('Raw (blue)','Adjusted (red)')
    end

    % Save output back into 4D file
    if is4D
        [pth,nm,ext] = fileparts(P(1,:));
        I = find(ext==',',1,'first');
        if ~isempty(I)
            ext = ext(1:I-1);
        end
        cspm_im_3Dto4D( char(Pd), [pth,filesep,'d',nm,ext] )
    end
end
% =========================================================================
function gt = cspm_globaltrend( P, pc, usespm, display )
    % CSPM_GLOBALTREND - Calculate and display global time trend of images.
    % T = CSPM_GLOBALTREND( P, PC, SPM, DISPLAY)
    %        P - cell array of files, e.g.,
    %            P = spm_get(Inf, '.img',{'Please select images for detrending'});
    %        PC - flag 1 = percent change relative to mean (default),
    %                  0 = absolute value
    %        SPM - flag 1 = use spm_global (all voxels above 1/8 of maximum
    %                       intensity; not very accurate)
    %                    0 = avearge of all voxels in volume (default)
    %        DISPLAY - flag 1 = make figure (default)
    %                       0 = no figure
    % Output: T - global time trend

    % @(#)cspm_globaltrend.m	1.0 Paul Macey 2005-08-01

    % Calculate global trend of series of images
    if nargin < 1
        if strcmp(spm('ver'),'SPM99') || strcmp(spm('ver'),'SPM2')
            P = spm_get(Inf, '.img','Select images for global trend calculation');
        else
            P = spm_select(Inf,'image','Select images for global trend calculation');
        end
        if isempty(P), return, end
    end
    if iscell(P)
        P = char(P);
    end
    if nargin < 2
        pc = 1;
    end
    if nargin < 3
        usespm = 0;
    end
    if nargin < 4
        display = 1;
    end

    V = spm_vol(P);
    nscan = length(V);
    gt = zeros(1,nscan);   % preallocate
    for i = 1:nscan
        if usespm
            % Note: spm_global estimates the mean after discounting voxels outside
            % the object using a criteria of greater than > (global mean)/8.
            % However, for large global signal changes, this does not accurately
            % estimate the global signal and LMGS detrending based on spm_global
            % does not remove all gobal components.
            gt(i) = spm_global(V(i));
        else
            Y = spm_read_vols(V(i));
            Y = Y(~isnan(Y));
            Y = Y(Y ~= 0);
            gt(i) = mean(mean(mean(Y)));
        end
    end
    if strcmp(spm('ver'),'SPM2')
        spm_close_vol(V);
    end

    if pc
        bl = mean(gt);
        gt = 100*(gt-bl)/bl;
        pstr = ' (% change)';
    else
        pstr = ' (raw)';
    end

    if display
        figure('Name',['Global timetrend for ',int2str(length(V)),' files',pstr])
        plot(gt)
        ylabel(['Signal Intensity',pstr])
        xlabel('Scan Number')
    end

    if nargout == 0
        clear gt
    end
end

% =========================================================================
function b = regression(y,X)
    % The following checks are omitted within LMGS (since we know what is being
    % passed to this function).
    %
    % Check that matrix (X) and left hand side (y) have compatible dimensions
    % [n,p] = size(X);
    % [n1,collhs] = size(y);
    % if n ~= n1,
    %     error('The number of rows in Y must equal the number of rows in X.');
    % end
    %
    % if collhs ~= 1,
    %     error('Y must be a vector, not a matrix');
    % end

    % Remove missing values, if any
    wasnan = (isnan(y) | any(isnan(X),2));
    if (any(wasnan))
        y(wasnan) = [];
        X(wasnan,:) = [];
    end

    % Find the least squares solution.
    [Q, R]=qr(X,0);
    b = R\(Q'*y);
end
% =========================================================================
function cspm_im_3Dto4D(P, outfile)
    % Convert series of 3D images to single 4D image volume (SPM5 only)
    % based on spm_config_3Dto4D, with some lint-informed changes.
    
    if nargin == 0
        P = spm_select(Inf,'image','Select images to group into 4D file');
        if isempty(P), return, end
        nfiles = size(P,1);
        Pnew = cell(1,nfiles);   % preallocate
        for i = 1:nfiles
            [pth,nm,ext] = fileparts(P(i,:));
            I = find(ext==',');
            if ~isempty(I)
                ext = ext(1:I-1);
            end
            Pnew{i} = [pth,filesep,nm,ext];
        end
        P = char(Pnew);
    elseif iscell(P)
        P = char(P);
    end
    if size(P,1) == 1
        msgbox('Only one image - selected multiple images for 4D file',...
            upper(mfilename),'warn')
        return
    end

    if nargin < 2
        [pth,nm,ext] = fileparts(P(1,:));
        [f,p] = uiputfile([pth,filesep,nm,'4D',ext],'Save as 4D file');
        if f == 0, return, end
        outfile = [p,f];
    end

    V    = spm_vol(P);
    ind  = cat(1,V.n);
    N    = cat(1,V.private);

    mx   = -Inf;
    mn   = Inf;
    for i=1:numel(V)
        dat      = V(i).private.dat(:,:,:,ind(i,1),ind(i,2));
        dat      = dat(isfinite(dat));
        mx       = max(mx,max(dat(:)));
        mn       = min(mn,min(dat(:)));
    end

    data_type = V(1).dt(1);
    r = spm_type(data_type,'maxval') - spm_type(data_type,'minval');
    if isfinite(r)
        sf = max(mx,-mn)/ r;
    else
        sf = 1;
    end
    ni         = nifti;
    ni.dat     = file_array(outfile,[V(1).dim numel(V)],data_type,0,sf,0);
    ni.mat     = N(1).mat;
    ni.mat0    = N(1).mat;
    ni.descrip = '4D image';
    create(ni);
    for i=1:size(ni.dat,4)
        ni.dat(:,:,:,i) = N(i).dat(:,:,:,ind(i,1),ind(i,2));
        spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
    end

    % Tidy up - delete 3D files
    for i=1:numel(V)
        delete(V(i).fname)
    end

end