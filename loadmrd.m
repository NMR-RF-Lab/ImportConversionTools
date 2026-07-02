function [dat] = loadmrd(MID, ~)
% LOADMRD - Import raw data *.dat files & header data. Requires FID-A software package.
%
% [data, info] = loadmrd(fID)
%
% Inputs
%         MID: number ID for measurement (must be in current folder)
%        keepos: *optional*; choice to keep os; 'keepos' or blank
%
% Outputs
%         dat: structure containing image & header

% Reformat MID to filename
MID = ['MID' num2str(MID) '.h5'];

% Check for file
if exist(MID, 'file')
    dset = ismrmrd.Dataset(MID, 'dataset');
else
    error(['File ' MID ' does not exist.  Please generate it.'])
end

% Load data
hdr = ismrmrd.xml.deserialize(dset.readxml);
D = dset.readAcquisition();

% Remove any noise scans
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
meas  = D.select(find(isNoise==0,1,'first'):D.getNumber);
clear D isNoise;

% Determine if 3D
if hdr.encoding.encodedSpace.matrixSize.z == 1
    dims = 2;
else
    dims = 3;
end

% Grab dimensional information
try nSli = hdr.encoding.encodingLimits.slice.maximum + 1; catch; nSli = 1; end
try nSet = hdr.encoding.encodingLimits.set.maximum + 1; catch; nSet = 1; end
try nSeg = hdr.encoding.encodingLimits.segment.maximum + 1; catch; nSeg = 1; end
try nEco = hdr.encoding.encodingLimits.contrast.maximum + 1; catch; nEco = 1; end
try nRep = hdr.encoding.encodingLimits.repetition.maximum + 1; catch; nRep = 1; end
try nAvg = hdr.encoding.encodingLimits.average.maximum + 1; catch; nAvg = 1; end
try nCha = hdr.acquisitionSystemInformation.receiverChannels; catch; nCha = 1; end
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Preallocate data
if dims == 3
    ksp_all = zeros([enc_Nx, enc_Ny, enc_Nz, nSli, nSet, nSeg, nEco, nRep, nAvg, nCha]);
    img = zeros([rec_Nx, rec_Ny, rec_Nz, nSli, nSet, nSeg, nEco, nRep, nAvg]);
    img_all = zeros([rec_Nx, rec_Ny, rec_Nz, nSli, nSet, nSeg, nEco, nRep, nAvg, nCha]);
else
    ksp_all = zeros([enc_Nx, enc_Ny, nSli, nSet, nSeg, nEco, nRep, nAvg, nCha]);
    img = zeros([rec_Nx, rec_Ny, nSli, nSet, nSeg, nEco, nRep, nAvg]);
    img_all = zeros([rec_Nx, rec_Ny, nSli, nSet, nSeg, nEco, nRep, nAvg, nCha]);
end

% Reconstruct images
for avg = 1:nAvg
    for rep = 1:nRep
        for eco = 1:nEco
            for seg = 1:nSeg
                for set = 1:nSet
                    for sli = 1:nSli
                        % Initialize the K-space storage array
                        K = zeros(enc_Nx, enc_Ny, enc_Nz, nCha);

                        % Select the appropriate measurements from the data
                        acqs = find((meas.head.idx.average == (rep-1)) ...
                            & (meas.head.idx.repetition == (rep-1)) ...
                            & (meas.head.idx.contrast == (eco-1)) ...
                            & (meas.head.idx.segment == (seg-1)) ...
                            & (meas.head.idx.set == (set-1)) ...
                            & (meas.head.idx.slice == (sli-1)));
                        for p = 1:length(acqs)
                            ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
                            kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1;
                            K(:,ky,kz,:) = meas.data{acqs(p)};
                        end

                        % Insert k-space into array
                        if dims == 3
                            ksp_all(:, :, :, sli, set, seg, eco, rep, avg, :) = K;
                        else
                            ksp_all(:, :, sli, set, seg, eco, rep, avg, :) = K;
                        end

                        % Reconstruct in x
                        K = fftshift(ifft(fftshift(K,1),[],1),1);

                        % Chop if needed
                        if (enc_Nx == rec_Nx)
                            im = K;
                        else
                            ind1 = floor((enc_Nx - rec_Nx)/2)+1;
                            ind2 = floor((enc_Nx - rec_Nx)/2)+rec_Nx;
                            im = K(ind1:ind2,:,:,:);
                        end

                        % Reconstruct in y then z
                        im = fftshift(ifft(fftshift(im,2),[],2),2);
                        if dims == 3
                            im = fftshift(ifft(fftshift(im,3),[],3),3);

                            % Insert into array
                            img_all(:, :, :, sli, set, seg, eco, rep, avg, :) = im;
                        else
                            img_all(:, :, sli, set, seg, eco, rep, avg, :) = im;
                        end

                        % Combine SOS across coils
                        im = sqrt(sum(abs(im).^2,4));

                        % Insert k-space into array
                        if dims == 3
                            img(:, :, :, sli, set, seg, eco, rep, avg) = im;
                        else
                            img(:, :, sli, set, seg, eco, rep, avg) = im;
                        end
                    end
                end
            end
        end
    end
end

% Combine into structure
dat.hdr = hdr;
dat.img = img;
dat.img_all = img_all;
dat.ksp_all = ksp_all;
end