% GIWAXS data analysis for semicrystalline conjugated polymer texture characterisation
% Hardware: Rigaku SmartLab (Lab 115, Clarendon Laboratory, Department of Physics, University of Oxford, UK)
% Version 1.0, by Bingjun Wang, 11/Oct/2019, Department of Physics, University of Oxford, UK
% Version 2.0, by Bingjun Wang, 25/Apr/2020, Department of Physics, University of Oxford, UK
% Version 2.1, by Bingjun Wang, 28/Apr/2020, Department of Physics, University of Oxford, UK
% Version 2.2, by Bingjun Wang, 06/Jun/2021, Department of Physics, University of Oxford, UK
% Email: bingjun.wang@physics.ox.ac.uk or bingjun.wang1995@outlook.com
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License: https://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Section 1: Variable Parameters - please change them for each sample

clearvars;
% About experiment data and settings
dataname = 'sample data';           % File name of the sample data, without .csv
bkgdataname = 'background data';	% File name of the substrate data, without .csv
alpha_c = 0.1500;                   % Sample critical angle, degree
alpha_i = 0.18;                     % Incident angle, degree

% About figure and data output
colourbarmin = 20;              % Minimum value of the colour bar in the exported figure, default = 20
colourbarmax = 800;             % Maximum value of the colour bar in the exported figure, default = 800
q_100 = 0.4;                    % Expected q value for (100) peak, angstrom^(-1)
delta_q_100 = 0.2;              % Width of the (100) q-cut arc, angstrom^(-1)
q_010 = 1.6;                    % Expected q-value for (010) linecut, angstrom^(-1)
delta_q_010 = 0.4;              % Width of the (010) q-cut arc, angstrom^(-1)
inplane_chi_min = 85;           % Minimum chi value of the in-plane cut, degree
inplane_chi_max = 90;           % Maximum chi value of the in-plane cut, degree
interval_100 = 0.5;             % Angle interval of exported (100) q-cut data, degree
interval_010 = 0.5;             % Angle interval of exported (010) q-cut data, degree
interval_inplane = 0.01;        % q interval of exported in-plane data, angstrom^(-1)
outputpath = 'file path';       % Full path for figure and data output


%% Section 2: Fixed Parameters - do not change them unless the diffractometer itself changes

DXL = 77.5;                     % Detector X length, mm
DYL = 38.5;                     % Detector Y length, mm
pixelsize = 0.1;                % Side length of a pixel, mm
D = 65;                         % Distance from detector centre to sample, mm
lambda = 1.540593;              % Cu K-alpha X-ray wavelength, angstrom


%% Section 3: Data Input and Preprocessing

M = csvread(strcat(dataname,'.csv'),2,0);
M_bkg = csvread(strcat(bkgdataname,'.csv'),2,0);
raw_X = M(:,1);                 % Pixel coordinate of X
raw_Y = M(:,2);                 % Pixel coordinate of Y
raw_I = M(:,3) - M_bkg(:,3);	% Intensity, background subtracted

% Make sure the all intensity values are meaningful in logarithmic scale
for i_1 = 1:1:size(M,1)
    if raw_I(i_1) < 0
        raw_I(i_1) = 0;
    end
end
I = raw_I + 1;                  % Intensity, background subtracted, ready to be shown in logarithmic scale


%% Section 4: Find Pixel Coordinates of the Horizontal Point M(X0, Y0)

% Find X0 - calculate collective X-intensity
X_intensity = zeros(DXL/pixelsize,2);
for i_2 = 0:1:(DXL/pixelsize - 1)
    X_intensity(i_2+1,1) = i_2;     % The first column: X-pixel index
    for i_3 = 0:1:(DYL/pixelsize - 1)
        X_intensity(i_2+1,2) = X_intensity(i_2+1,2) + M(i_2+i_3*DXL/pixelsize+1,3);     % The second column: total intensity for each X-pixel
    end
end

% Find X0 - find all local maxima of intensity
X_localmax = zeros(DXL/pixelsize,2);
for i_4 = 2:1:(DXL/pixelsize-1)
    if X_intensity(i_4,2) > X_intensity(i_4+1,2) && X_intensity(i_4,2) > X_intensity(i_4-1,2)
        X_localmax(i_4,1) = X_intensity(i_4,1);
        X_localmax(i_4,2) = X_intensity(i_4,2);
    end
end

% Find X0 - delete blank rows, sort rows, and finally find X0 (in pixel coordinate)
X_localmax(all(X_localmax == 0,2),:) = [];
X_localmax = sortrows(X_localmax,-2);
X0_pixel = 0.5*(X_localmax(1,1) + X_localmax(2,1));

% Find Y0 - calculate collective Y-intensity
safe_Y = 20;        % Set range for Y0; the pixel coordinate of Y0 will be absolutely smaller than this value
Y_intensity = zeros(safe_Y+1,2);
for i_5 = 0:1:safe_Y
    Y_intensity(i_5+1,1) = i_5;     % The first column: Y-pixel index
    for i_6 = 0:1:(0.4*DXL/pixelsize - 1)         % Exclude the incluence close to the specular reflection point
        Y_intensity(i_5+1,2) = Y_intensity(i_5+1,2) + M(i_5*DXL/pixelsize+i_6+1,3);     % The second column: total intensity for each Y-pixel
    end
end

% Find Y0 - sort rows, find Yv (in pixel coordinate), calaulate Y0 (in pixel coordinate) and gamma
Y_intensity = sortrows(Y_intensity,-2);
Yv_pixel = Y_intensity(1,1);
Y0_pixel = (D*Yv_pixel*pixelsize - (D*D + 0.25*DYL*DYL - 0.5*DYL*Yv_pixel*pixelsize) * tan(deg2rad(alpha_c)))/(D*pixelsize + (Yv_pixel*pixelsize - 0.5*DYL)*pixelsize*tan(deg2rad(alpha_c)));
gamma = atan((0.5*DYL - Y0_pixel*pixelsize) / D);       % Detector rotation angle, rad


%% Section 5: Conversion from Pixel Coordinate to 3-D Coordinate

td = zeros(size(M,1),3);        % The three columns represent x, y, and z coordinates in 3-D coordinate, respectively
td(:,1) = abs((M(:,1) - X0_pixel) * pixelsize);                                 % 3-D x coordinates
td(:,2) = D/cos(gamma) - (M(:,2) - Y0_pixel) * pixelsize * sin(gamma);          % 3-D y coordinates
td(:,3) = (M(:,2) - Y0_pixel) * pixelsize * cos(gamma);                         % 3-D z coordinates


%% Section 6: Extraction of alpha_f and two_theta_f

angles = zeros(size(M,1),2);            % The first column is alpha_f (rad), and the second column is two_theta_f (rad)

for i_7 = 1:1:size(M,1)
    angles(i_7,1) = pi/2 - acos(td(i_7,3) / sqrt(td(i_7,1)*td(i_7,1) + td(i_7,2)*td(i_7,2) + td(i_7,3)*td(i_7,3)));        % alpha_f, rad
    angles(i_7,2) = acos(td(i_7,2) / sqrt(td(i_7,1)*td(i_7,1) + td(i_7,2)*td(i_7,2)));          %two_theta_f, rad
end


%% Section 7: Calculation of Components of Scattering Vector (q_x, q_y, q_xy, and q_z)

q = zeros(size(M,1),4);
q(:,1) = 2*pi/lambda * cos(angles(:,1)) .* sin(angles(:,2));                                 % q_x, A^(-1)
q(:,2) = 2*pi/lambda * (cos(angles(:,1)) .* cos(angles(:,2)) - cos(deg2rad(alpha_i)));       % q_y, A^(-1)
q(:,3) = sqrt(q(:,1) .* q(:,1) + q(:,2) .* q(:,2));                                          % q_xy, A^(-1)
q(:,4) = 2*pi/lambda * (sin(angles(:,1)) + sin(deg2rad(alpha_i)));                           % q_z, A^(-1)


%% Section 8: q_z - q_xy Plot

% Correction for wrong peak near the centre
for i_8=1:1:size(M,1)
    if q(i_8,3) < 0.12 && q(i_8,4) < 0.25
        I(i_8) = 1;                          % Ensure the intensity values are meaningful in log scale
    end
end

% Plot the figure
set(0,'defaultfigurecolor','w');
tri=delaunay(q(:,3),q(:,4));
patch('Faces',tri,'Vertices',[q(:,3),q(:,4)],'FaceVertexCData',I,'FaceColor','interp','EdgeColor','none');
caxis([colourbarmin,colourbarmax]);
colorbar;
colormap jet;
set(gca,'ColorScale','log');
backColor = [0 0 0];
set(gca,'color',backColor);
axis equal;
axis([0 2 0 2]);
set(gca,'xtick',0:0.4:2,'ytick',0:0.4:2,'tickdir','out','fontname','cmr12','fontsize',12,'fontweight','bold');
xlabel('\boldmath $q_{xy} \ (\mathrm{\AA}^{-1})$','interpreter','latex','fontsize',14);
ylabel('\boldmath $q_{z} \ (\mathrm{\AA}^{-1})$','interpreter','latex','fontsize',14);
set(gca,'units','normalized','position',[0.15,0.15,0.8,0.8]);

figure(1);
set(gcf,'InvertHardCopy','off');


%% Section 9: q-cut Data Production (for further analysis in Origin)

% Matrix definition
M_100_raw = zeros(size(M,1),2);             % First column: chi (degree); second column: intensity
M_100 = zeros(90/interval_100 + 1,4);       % First column: chi (degree); second column: intensity; third column: No. of points; fourth column: average intensity
M_010_raw = zeros(size(M,1),2);             % First column: chi (degree); second column: intensity
M_010 = zeros(90/interval_010 + 1,4);       % First column: chi (degree); second column: intensity; third column: No. of points; fourth column: average intensity
q_chi = zeros(size(M,1),2);                 % First column: q (A^-1); second column: chi (degree)

% Parameter definition
q_100_min = q_100 - 0.5 * delta_q_100;
q_100_max = q_100 + 0.5 * delta_q_100;
q_010_min = q_010 - 0.5 * delta_q_010;
q_010_max = q_010 + 0.5 * delta_q_010;

% Calculation of q and chi
q_chi(:,1) = sqrt(q(:,1) .* q(:,1) + q(:,2) .* q(:,2) + q(:,4) .* q(:,4));      % q, A^(-1)
q_chi(:,2) = rad2deg(atan(q(:,3)./q(:,4)));                                     % chi, degree

for i_9 = 1:1:(90/interval_100 + 1)                                             % chi values for output
    M_100(i_9,1) = (i_9 - 1) * interval_100;
    M_010(i_9,1) = M_100(i_9,1);
end

% Intensity accumulation
i_10 = 1;
for i_11 = 1:1:size(M,1)
    if q_chi(i_11,1) > q_100_min && q_chi(i_11,1) <= q_100_max && q_chi(i_11,2) > 0
        M_100_raw(i_10,1) = q_chi(i_11,2);
        M_100_raw(i_10,2) = raw_I(i_11);
        i_10 = i_10 + 1;
    end
end
M_100_raw(all(M_100_raw == 0,2),:) = [];        % Delete all blank rows
M_100_raw = sortrows(M_100_raw);                % Sort the matrix in ascending order of chi

i_12 = 1;
for i_13 = 1:1:size(M,1)
    if q_chi(i_13,1) > q_010_min && q_chi(i_13,1) <= q_010_max && q_chi(i_13,2) > 0
        M_010_raw(i_12,1) = q_chi(i_13,2);
        M_010_raw(i_12,2) = raw_I(i_13);
        i_12 = i_12 + 1;
    end
end
M_010_raw(all(M_010_raw == 0,2),:) = [];        % Delete all blank rows
M_010_raw = sortrows(M_010_raw);                % Sort the matrix in ascending order of chi

% Merge data into every chi interval
i_13 = 1;
for i_14 = 1:1:(90/interval_100 + 1)
    for i_13 = i_13:1:size(M_100_raw,1)
        if M_100_raw(i_13,1) > M_100(i_14,1) - 0.5*interval_100 && M_100_raw(i_13,1) <= M_100(i_14,1) + 0.5*interval_100
            M_100(i_14,2) = M_100(i_14,2) + M_100_raw(i_13,2);      % Intensity accumulation of each angle interval
            M_100(i_14,3) = M_100(i_14,3) + 1;                      % Counter of No. of points
        else
            break
        end
    end
end
M_100(:,4) = M_100(:,2) ./ M_100(:,3);          % Calculate average intensity
M_100(any(M_100(:,3) == 0,2),:) = [];           % Delete all rows that have no data (i.e. No. of points = 0)

i_15 = 1;
for i_16 = 1:1:(90/interval_010 + 1)
    for i_15 = i_15:1:size(M_010_raw,1)
        if M_010_raw(i_15,1) > M_010(i_16,1) - 0.5*interval_010 && M_010_raw(i_15,1) <= M_010(i_16,1) + 0.5*interval_010
            M_010(i_16,2) = M_010(i_16,2) + M_010_raw(i_15,2);      % Intensity accumulation of each angle interval
            M_010(i_16,3) = M_010(i_16,3) + 1;                      % Counter of No. of points
        else
            break
        end
    end
end
M_010(:,4) = M_010(:,2) ./ M_010(:,3);          % Calculate average intensity
M_010(any(M_010(:,3) == 0,2),:) = [];           % Delete all rows that have no data (i.e. No. of points = 0)


%% Section 10: In-Plane Data Production (for further analysis in Origin)

% Matrix definition
M_inplane_raw = zeros(size(M,1),2);                 % First column: q (A^(-1)); second column: intensity
M_inplane = zeros(2.5/interval_inplane + 1,4);      % First column: q (A^(-1)); second column: intensity; third column: No. of points; fourth column: average intensity

for i_17 = 1:1:(2.5/interval_inplane + 1)           % q values for output
    M_inplane(i_17,1) = (i_17 - 1) * interval_inplane;
end

% Intensity accumulation
i_18 = 1;
for i_19 = 1:1:size(M,1)
    if abs(q_chi(i_19,2)) > inplane_chi_min && abs(q_chi(i_19,2)) <= inplane_chi_max
        M_inplane_raw(i_18,1) = q_chi(i_19,1);
        M_inplane_raw(i_18,2) = raw_I(i_19);
        i_18 = i_18 + 1;
    end
end
M_inplane_raw(all(M_inplane_raw == 0,2),:) = [];        % Delete all blank rows
M_inplane_raw = sortrows(M_inplane_raw);                % Sort the matrix in ascending order of q

% Merge data into every q interval
i_20 = 1;
for i_21 = 1:1:(2.5/interval_inplane + 1)
    for i_20 = i_20:1:size(M_inplane_raw,1)
        if M_inplane_raw(i_20,1) > M_inplane(i_21,1) - 0.5*interval_inplane && M_inplane_raw(i_20,1) <= M_inplane(i_21,1) + 0.5*interval_inplane
            M_inplane(i_21,2) = M_inplane(i_21,2) + M_inplane_raw(i_20,2);          % Intensity accumulation of each angle interval
            M_inplane(i_21,3) = M_inplane(i_21,3) + 1;                              % Counter of No. of points
        else
            break
        end
    end
end
M_inplane(:,4) = M_inplane(:,2) ./ M_inplane(:,3);          % Calculate average intensity
M_inplane(any(M_inplane(:,3) == 0,2),:) = [];               % Delete all rows that have no data (i.e. No. of points = 0)


%% Section 11: Figure and Data Output

% Filename definition
q_100_name = strcat(dataname,'-q cut-',num2str(q_100),'-',num2str(delta_q_100));       % Filename of the exported csv file of (100) q-cut
q_010_name = strcat(dataname,'-q cut-',num2str(q_010),'-',num2str(delta_q_010));       % Filename of the exported csv file of (010) q-cut
figure_name = strcat(dataname,'-',num2str(colourbarmin),'-',num2str(colourbarmax));    % Filename of the exported figure

% Convert matrices to tables
t_M_100 = array2table(M_100,'VariableNames',{'chi','Intensity','No. of Points','Average Intensity'});
t_M_010 = array2table(M_010,'VariableNames',{'chi','Intensity','No. of Points','Average Intensity'});
t_M_inplane = array2table(M_inplane,'VariableNames',{'q','Intensity','No. of Points','Average Intensity'});

% Output
filepath = pwd;
cd(outputpath);
print(gcf,strcat(figure_name,'.tif'),'-dtiffnocompression', '-r300');
writetable(t_M_100,strcat(q_100_name,'.csv'));
writetable(t_M_010,strcat(q_010_name,'.csv'));
writetable(t_M_inplane,strcat(dataname,'-in-plane cut.csv'));
cd(filepath);
