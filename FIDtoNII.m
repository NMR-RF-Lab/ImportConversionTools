function nii = FIDtoNII(fname, img, p_full)
% Converts an imported *.fid file to a *.nii file.

% Inputs
%    fname: file name
%      img: complex image data
%   p_full: parameters obtained from procpar file

% Outputs
%      nii: NIfTI structure containing both image & header data

% Specify origin
vs = [(p_full.lro*10)/(p_full.np/2) (p_full.lpe*10)/(p_full.nv) (p_full.thk)]; % voxel size

% Create description
if length(fname) == 2
    des = char(regexprep(fname(1), '_', ' '));
    des = ['grmd ', des(1:end-7), des(end-3:end), '2022'];
else
    des = ['grmd ', regexprep(fname, '_', ' '), '2022']; 
    % Remove 'p1'
    
end

% Generate nii file in matlab
nii = make_nii(img, vs, [], [], des);

% Make changes to header
nii.hdr.dime.xyzt_units = 2; % change units to millimeter

% Save file
scheck = input('Would you like to save the file? [Y/N] ', 's');
if strcmp(scheck, 'Y')
    fname_nii = input('Enter the file name (with file extension): ', 's');
    if ~strcmp(fname_nii(end-3), '.')
        error('A file extension was not included in the file name.')
    end
    comp_check = input('Do you want the data saved in complex form? [Y/N] ', 's');
    switch comp_check
        case 'Y'
            nii.hdr.dime.datatype = 1792; nii.hdr.dime.bitpix = 1792; % change datatype to complex double
            save_nii(nii, fname_nii);

        case 'N'
            fname_nii_mag = [fname_nii(1:end-4), '_Mag', fname_nii(end-3:end)];
            fname_nii_pha = [fname_nii(1:end-4), '_Phase', fname_nii(end-3:end)];

            nii_mag = nii;
            nii_mag.img = abs(nii.img);
            nii_mag.hdr.dime.datatype = 64; nii_mag.hdr.dime.bitpix = 64;
            save_nii(nii_mag, fname_nii_mag);

            nii_pha = nii;
            nii_pha.img = angle(nii.img);
            nii_pha.hdr.dime.datatype = 64; nii_pha.hdr.dime.bitpix = 64;
            save_nii(nii_pha, fname_nii_pha);
    end
end