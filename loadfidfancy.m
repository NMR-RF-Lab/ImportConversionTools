function [K, I, P, p, F] = loadfidfancy(varargin)
%LOADFIDFANCY Imports 2D & 3D *.fid folders & corresponding parameters. Requires "readprocpar & "loadfid" or "loadfid2" to be on path.
%
%   [K, I, P, p, F] = loadfidfancy(fname, Narea, cshift, zip)
%
%   Inputs (name-value arguments)
%             fid : *optional* *.fid folder name & path. If not inputted, function will open dialog box & request user input
%       DCremoval : *optional* removes DC artifact, based on given input area; default = 0.2
%   circularshift : *optional* circular shift peak k-space value to center; default = false
%             zip : *optional* 2x zero-interpolation filling multiplier; default = false
%  
%   Outputs
%               K : k-space matrix
%               I : image matrix
%               P : important metadata parameters
%               p : all metadata parameters in "readprocpar" format
%               F : file name & path
%
% Jacob Degitz, Texas A&M University
% Created 7/13/2022
% Last edited 4/11/2024

% Check for each input
CHECK_fid = true;
CHECK_DCr = true;
CHECK_csh = true;
CHECK_zip = true;
for i = 1:nargin
    % Check for folder input
    if strcmp(varargin{i}, 'fid') && CHECK_fid
        % Remove input from options
        CHECK_fid = false;

        % Ensure it is a string or character
        if ischar(varargin{i+1}) || isstring(varargin{i+1})
            % Save variable
            fid = char(varargin{i+1});
        else
            error('Folder path must be either a string or character array.')
        end

    % Check for DCremoval input
    elseif strcmp(varargin{i}, 'DCremoval') && CHECK_DCr
        % Remove input from options
        CHECK_DCr = false;

        % Ensure it is a numeric > 0 & < 1
        if ~isnumeric(varargin{i+1})
            error('DC removal area must be number.')
        elseif varargin{i+1} <= 0
            error('DC removal area must greater than zero.')
        elseif varargin{i+1} >= 1
            error('DC removal area must less than one.')
        else
            % Save variable
            Narea = varargin{i+1};
        end

    % Check for circularshift input
    elseif strcmp(varargin{i}, 'circularshift') && CHECK_csh
        % Remove input from options
        CHECK_csh = false;

        % Ensure it is a logical
        if ~islogical(varargin{i+1})
            error('circularshift option must be logical.')
        else
            % Save variable
            cshiftFLAG = varargin{i+1};
        end

    % Check for zip input
    elseif strcmp(varargin{i}, 'zip') && CHECK_zip
        % Remove input from options
        CHECK_zip = false;

        % Ensure it is a logical
        if ~islogical(varargin{i+1})
            error('zip option must be logical.')
        else
            % Save variable
            zipFLAG = varargin{i+1};
        end
    end   

    % Shift forward indexer
    i = i+1;
end

% Add default options based on which inputs weren't grabbed
if CHECK_fid
    fid = uigetdir();
end
if CHECK_DCr
    Narea = 0.2;
end
if CHECK_csh
    cshiftFLAG = false;
end
if CHECK_zip
    zipFLAG = false;
end

% Obtain file path 
if strcmp(fid(end-3:end), '.fid')
    F = what(fid);
else
    F = what([fid '.fid']);
end

% Isolate file name from path
[F, fid, ~] = fileparts(F.path);

% Combine if not in current directory
if ~strcmp(pwd, F)
    fid = [F '\' fid];
end

% Load data, parameters
% Check if data is p1 of two multislice scans
if strcmp(fid(end-6:end-5), 'p1') == 1
    %%% Combine files
    fname1 = fid;

    % Obtain second file name
    fname2 = fname1(end-6:end); % "p1" & date
    fname2 = split(fname2, "_"); % split "p1" from date
    fname2 = char(strcat(fname1(1:end-7), "p2_", string(fname2(2)))); % recombine date with "p2" & filename to create p2 filename

    % Load Data & Parameters
    K1 = loadfid_multislice(fname1);
    K2 = loadfid_multislice(fname2);

    %%% Import both sets, everything should be the same except for slice offset
    p = readprocpar([fname1, '.fid']);
    p2 = readprocpar([fname2, '.fid']);
    p.ns = p.ns*2; % double the number of slices
    p.pss = [p.pss p2.pss]; % add slice offsets from second file
else
    % Second check for multislice (for only 1 file)
    if contains(fid, 'mult')
        kspace_unf = loadfid_multislice(fid);
        p = readprocpar([fid, '.fid']);
    else
        % Check for either loadfid or loadfid2
        if exist('loadfid', 'file') == 2
            kspace_unf = loadfid(fid); % Load k-space data & parameters
        elseif exist('loadfid2', 'file') == 2
            kspace_unf = loadfid2(fid);
        end
        p = readprocpar([fid, '.fid']);

        % Check if 3D image not using 'mult' in title
        if p.ns > 1
            clear kspace_unf
            kspace_unf = loadfid_multislice(fid);
        end
    end
end

% Not sure if 'dimX/Y/Z' or 'posX/Y/Z' is location of gradients
if isequal(p.dimX, {'lpe2'})
    p.dimX = {'gss'};
    if isequal(p.dimY, {'lpe'})
        p.dimY = {'gpe'}; p.dimZ = {'gro'};
    else
        p.dimY = {'gro'}; p.dimZ = {'gpe'};
    end
elseif isequal(p.dimY, {'lpe2'})
    p.dimY = {'gss'};
    if isequal(p.dimX, {'lpe'})
        p.dimX = {'gpe'}; p.dimZ = {'gro'};
    else
        p.dimX = {'gro'}; p.dimZ = {'gpe'};
    end
elseif isequal(p.dimZ, {'lpe2'})
    p.dimZ = {'gss'};
    if isequal(p.dimX, {'lpe'})
        p.dimX = {'gpe'}; p.dimY = {'gro'};
    else
        p.dimX = {'gro'}; p.dimY = {'gpe'};
    end
end

% Isolate important parameters
P.Info = struct('date',p.date,'seq',p.seqfil,'B0',p.B0/10000,'gax',struct('gX',p.dimX,'gY',p.dimY,'gZ',p.dimZ)); % date, pulse sequence, magnet field strength (T), gradient along x-, y- & z-axis
P.ImagParam = struct('tr',p.tr,'te',p.te,'tped',p.tped,'nechos',p.ne); % repetition time (s), echo time (s), phase encode time (s), # echoes
P.Matrix = struct('nro',p.np/2,'npe',p.nv); % # readout points, # phase points
P.FOV = struct('fov',struct('ro',p.lro,'pe',p.lpe),'offro',p.pro,'orient',p.orient); % readout FOV (cm), phase FOV (cm), readout offset (cm), orientation
P.SliceInfo = struct('thks',p.thk/10,'offs',p.pss,'ns',p.ns); % slice thickness (cm), slice offset (cm), # slices
P.Coils = struct('coilg',p.gcoil,'coilrf',p.rfcoil); % gradient coil, rf coil
P.RF = struct('fa',p.fliplist,'dur',struct('p90',p.p1,'p180',p.p2),'ptn',struct('p90',p.p1pat,'p180',p.p2pat),'pwr',struct('p90',p.tpwr1,'p180',p.tpwr2)); % pulse duration (us), pattern, power (dB)
P.Acquisition = struct('freq',p.sfrq,'sw',p.sw,'tacq',p.at,'recgain',p.gain); % observe frequency (MHz), spectral width (Hz), acquisition time (ms), receiver gain (dB)
P.Grad = struct('gro',p.gro,'gpe',p.gpe,'gss',p.gss,'gssr',p.gssr);

% Combining if multislice
if strcmp(fid(end-6:end-5), 'p1') == 1
    kspace_unf = zeros(P.Matrix.nro, P.Matrix.npe, p.ns);

    for j = 1:p.ns
        % Determine if first or second file
        if rem(j,2) == 1 
            kspace_unf(:,:,j) = K1(:,:,1); % add slice to main matrix
            % Check if iteration is not last
            if j < p.ns-1
                K1 = K1(:,:,2:end); % subtract slice from individual matrix
            end
        else
            kspace_unf(:,:,j) = K2(:,:,1); % add slice to main matrix
            % Check if iteration is not last
            if j < p.ns
                K2 = K2(:,:,2:end); % subtract slice from individual matrix
            end
        end
    end
end

% Check if k-space data will be zipped
if zipFLAG
    % Update parameters
    zfact = 2;
    nro = P.Matrix.nro;
    npe = P.Matrix.npe;
    P.Matrix.nro = P.Matrix.nro*zfact;
    P.Matrix.npe = P.Matrix.npe*zfact;
end

% Preallocate data
I = zeros(P.Matrix.nro, P.Matrix.npe, P.SliceInfo.ns);
K = zeros(P.Matrix.nro, P.Matrix.npe, P.SliceInfo.ns);

% Basic filtering
for j = 1:P.SliceInfo.ns
    % Filter in frequency domain
    kspace_temp = filt_FD(kspace_unf(:,:,j), size(kspace_unf, 1), size(kspace_unf, 2), Narea);

    % Circular shift
    if cshiftFLAG
        kspace_temp = cent(kspace_temp);
    end

    % ZIP
    if zipFLAG
        kspace_temp = padarray(kspace_temp, [(P.Matrix.nro-nro)/2, (P.Matrix.npe-npe)/2], 0, 'both');
    end

    % Filter in time domain
    img_temp = filt_TD(kspace_temp, size(kspace_temp, 1), size(kspace_temp, 2), floor(Narea/2));

    kspace_temp = ifftshift(fft2(ifftshift(img_temp)));

    K(:,:,j) = kspace_temp;
    I(:,:,j) = img_temp;
end

    function img = filt_TD(kspace, nro, npe, N_factor)
        % Basic time domain filtering 

        % Area for obtaining mean kspace background
        cols = floor(N_factor*npe);

        img = fftshift(ifft2(fftshift(kspace)));

        % Remove noise line in true image
        rng(90); randsign1 = (rand(nro,1) > 0.5)*2 - 1; rng(12); randsign2 = (rand(nro,1) > 0.5)*2 - 1;
        [datastd_r, datamn_r] = std(real(img(:,2:2+floor(cols/2))), 0, 2); % noise std & mean (real)
        [datastd_i, datamn_i] = std(imag(img(:,2:2+floor(cols/2))), 0, 2); % noise std & mean  (imag)
        rand_real = datamn_r + datastd_r.*randsign1; % approximated first column (real)
        rand_imag = datamn_i + datastd_i.*randsign2; % approximated first column (imag)
        img(:,1) = complex(rand_real, rand_imag);

        % Check if DC noise is still present (if so, replace w/ smaller value)
        dif = 1.3;
        DC = [(size(img,1)/2)+1, (size(img,2)/2)+1];
        if abs(img(DC(1),DC(2))) > (mean2(abs(img(DC(1)-1:DC(1)+1,DC(2)-1:DC(2)+1)))*dif)
            img(DC(1),DC(2)) = abs(img(DC(1),DC(2))/dif)*exp((1i*angle(img(DC(1),DC(2)))));
        end
    end

    function kspace_cs = cent(kspace)
        % Circular shifts k-space data (centers k-space & helps with unwrapping)

        dims = size(kspace);

        % Check if 2D
        if size(dims, 2) == 2
            dims(1,3) = 1;
            [pk_row, pk_col] = find(abs(kspace(:,:))==max(max(abs(kspace(:,:))))); % Identify offset
        else
            [pk_row, pk_col] = find(abs(kspace(:,:,1))==max(max(max(abs(kspace(:,:,1))))));
        end

        kspace_cs = zeros(size(kspace(:,:,:)));

        for z = 1:dims(3) % Circular shift all kspace data
            kspace_cs(:,:,z) = circshift(kspace(:,:,z), (dims(1)/2) - pk_row, 1);
            kspace_cs(:,:,z) = circshift(kspace_cs(:,:,z), (dims(2)/2) - pk_col, 2);
        end
    end

    function kspace = filt_FD(kspace_unf, nro, npe, N_factor)
        % Basic Fourier domain filtering

        % Area for obtaining mean kspace background
        rows = floor(N_factor*nro);
        cols = floor(N_factor*npe);
        datamn = mean2(kspace_unf(end-rows:end, end-cols:end)); % mean noise value
        kspace_unf = kspace_unf - datamn; % subtract mean noise from image

        % Remove noise line in k-space
        kspace = kspace_unf;
        if max(abs(kspace_unf(:,1))) > mean2(abs(kspace_unf))
            rng(88); randsign1 = (rand(nro,1) > 0.5)*2 - 1; rng(74); randsign2 = (rand(nro,1) > 0.5)*2 - 1;
            [datastd_r, datamn_r] = std(real(kspace_unf(:,1:1+floor(cols/2))), 0, 2); % noise std & mean (real)
            [datastd_i, datamn_i] = std(imag(kspace_unf(:,1:1+floor(cols/2))), 0, 2); % noise std & mean  (imag)
            rand_real = datamn_r + datastd_r.*randsign1; % approximated first column (real)
            rand_imag = datamn_i + datastd_i.*randsign2; % approximated first column (imag)
            kspace(:,1) = complex(rand_real, rand_imag);
        end
    end

    function kspace_unfilt = loadfid_multislice(imagefile)
        % This function seperates one kspace matrix into the number of slices you obtained
        % Order should be changed if there are more slices

        % Checks for readprocpar function
        if exist('readprocpar', 'file') == 2
            prms = readprocpar(sprintf('%s.fid',imagefile));
        else
            error('"readprocpar" cannot be found.')
        end

        [~,order] = sort(prms.pss,'ascend'); % sort data by ascending order

        % Checks for loadfid or loadfid2
        if exist('loadfid', 'file') == 2
            kspace_multislice = loadfid(imagefile);
        elseif exist('loadfid2', 'file') == 2
            kspace_multislice = loadfid2(imagefile);
        else
            error('Both "loadfid" and "loadfid2" cannot be found.')
        end

        n = length(order);

        % Initialize empty arrays
        kspace_unfilt = zeros(prms.np/2, prms.nv, n);

        for z = 1:n
            % Load all slices
            kspace_unfilt(:,:,z) = kspace_multislice(:,order(z):n:end);
        end
    end
end