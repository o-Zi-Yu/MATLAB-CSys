fc = 25e9; % carrier frequency (Hz)
bsPosition = [42.729727, -73.680000]; % lat, lon
bsAntSize = [8 8]; % number of rows and columns in rectangular array (base station)
bsArrayOrientation = [-30 0].'; % azimuth (0 deg is East, 90 deg is North) and elevation (positive points upwards) in deg
ueAntSize = [2 2]; % number of rows and columns in rectangular array (UE).

uePosition = [42.728698, -73.680831]; % lat, lon
ueArrayOrientation = [180 45].';

uePosition2 = [42.729165, -73.681293];
ueArrayOrientation2 = [90 30].'; % azimuth (0 deg is East, 90 deg is North) and elevation (positive points upwards)  in deg

uePosition3 = [42.728801, -73.679098];
ueArrayOrientation3 = [45 60].'; % azimuth (0 deg is East, 90 deg is North) and elevation (positive points upwards)  in deg

uePosition4 = [42.729503, -73.679544];
ueArrayOrientation4 = [0 0].'; % azimuth (0 deg is East, 90 deg is North) and elevation (positive points upwards)  in deg

reflectionsOrder = 6; % number of reflections for ray tracing analysis (0 for LOS)

% Bandwidth configuration, required to set the channel sampling rate and for perfect channel estimation
SCS = 30; % subcarrier spacing
NRB = 52; % number of resource blocks, 10 MHz bandwidth
nSlot = 0;
nLayers = 1;
scOffset = 0; % no offset
noRBs = 1; % average channel conditions over 1 RB to calculate beamforming weights

% Set carrier resource grid properties (30 kHz SCS and 20 MHz bandwidth)
%carrier = nrCarrierConfig;
carrier.NCellID = 0;
carrier.SubcarrierSpacing = 30;
carrier.CyclicPrefix = "normal";
carrier.NSizeGrid = 52;
carrier.NStartGrid = 0;
carrier = redcap_CarrierConfig(carrier);

% Set number of transmit and receive antennas
simParameters = struct; % Create simParameters structure
simParameters.NFrames = 5; % Number of 10 ms frames
simParameters.SNRIn = 15; % SNR range (dB)

simParameters.NTxAnts = 1;
simParameters.NRxAnts = 1;

if exist('viewer', 'var') && isvalid(viewer) % viewer handle exists and viewer window is open
    viewer.clearMap();
else
    viewer = siteviewer("Basemap", "openstreetmap", "Buildings", "map.osm");
end

%create bs and ue location

bsSite = txsite("Name", "Base station", ...
"Latitude", bsPosition(1), "Longitude", bsPosition(2), ...
    "AntennaAngle", bsArrayOrientation(1:2), ...
    "AntennaHeight", 4, ... % in m
    "TransmitterFrequency", fc);
%bsSite.Antenna.Taper = 0;
ueSite = rxsite("Name", "UE", ...
"Latitude", uePosition(1), "Longitude", uePosition(2), ...
    "AntennaHeight", 1, ... % in m
    "AntennaAngle", ueArrayOrientation(1:2));

%visulize the locations

bsSite.show();
ueSite.show();

%analysis of signal reciving range without reflection
rtpm = propagationModel("raytracing", "Method", "sbr", MaxNumReflections = 0, ...
BuildingsMaterial = "perfect-reflector", TerrainMaterial = "perfect-reflector");
coverage(bsSite, rtpm, "SignalStrengths", -80:-5)

%ray tracing analysis
pm = propagationModel("raytracing", "Method", "sbr", "MaxNumReflections", reflectionsOrder);

waveformInfo = redcap_OFDMInfo(carrier.NSizeGrid, carrier.SubcarrierSpacing, carrier.CyclicPrefix);

c = physconst('LightSpeed');
lambda = c / fc;
ueArray = phased.NRRectangularPanelArray('Size', [ueAntSize(1:2) 1 1], 'Spacing', [0.5 * lambda * [1 1] 1 1]);
ueArray.ElementSet = {phased.IsotropicAntennaElement}; % isotropic antenna element
bsArray = phased.NRRectangularPanelArray('Size', [bsAntSize(1:2) 1 1], 'Spacing', [0.5 * lambda * [1 1] 1 1]);
bsArray.ElementSet = {phased.NRAntennaElement('PolarizationAngle', -45) phased.NRAntennaElement('PolarizationAngle', 45)}; % cross polarized elements

rays = raytrace(bsSite, ueSite, pm, "Type", "pathloss");

try
    plot(rays{1})

    pathToAs = [rays{1}.PropagationDelay] - min([rays{1}.PropagationDelay]); % Time of arrival of each ray (normalized to 0 sec)
    avgPathGains =- [rays{1}.PathLoss]; % Average path gains of each ray
    pathAoDs = [rays{1}.AngleOfDeparture]; % AoD of each ray
    pathAoAs = [rays{1}.AngleOfArrival]; % AoA of each ray
    isLOS = any([rays{1}.LineOfSight]);

    channel = nrCDLChannel;
    channel.DelayProfile = 'Custom';
    channel.PathDelays = pathToAs;
    channel.AveragePathGains = avgPathGains;
    channel.AnglesAoD = pathAoDs(1, :); % azimuth of departure
    channel.AnglesZoD = 90 - pathAoDs(2, :); % channel uses zenith angle, rays use elevation
    channel.AnglesAoA = pathAoAs(1, :); % azimuth of arrival
    channel.AnglesZoA = 90 - pathAoAs(2, :); % channel uses zenith angle, rays use elevation
    channel.HasLOSCluster = isLOS;
    channel.CarrierFrequency = fc;
    channel.NormalizeChannelOutputs = false; % do not normalize by the number of receive antennas, this would change the receive power
    channel.NormalizePathGains = false; % set to false to retain the path gains

    channel.SampleRate = waveformInfo.SampleRate;

    % UE array (single panel)
    channel.ReceiveAntennaArray = ueArray;
    channel.ReceiveArrayOrientation = [ueArrayOrientation(1); (-1) * ueArrayOrientation(2); 0]; % the (-1) converts elevation to downtilt
    % Base station array (single panel)
    channel.TransmitAntennaArray = bsArray;
    channel.TransmitArrayOrientation = [bsArrayOrientation(1); (-1) * bsArrayOrientation(2); 0]; % the (-1) converts elevation to downtilt
    channel.ChannelFiltering = false;
    [pathGains, sampleTimes] = channel();

    pg = permute(pathGains, [2 1 3 4]); % first dimension is the number of paths

    if isLOS
        % in LOS cases sum the first to paths, they correspond to the LOS ray
        pg = [sum(pg(1:2, :, :, :)); pg(3:end, :, :, :)];
    end

    pg = abs(pg).^2;

    % plot(pow2db(pg(:, 1, 1, 1)), 'o-.'); hold on
    % plot(avgPathGains, 'x-.'); hold off
    % legend("Instantaneous (1^{st} tx - 1^{st} rx antenna)", "Average (from ray tracing)")
    % xlabel("Path number"); ylabel("Gain (dB)")
    % title('Path gains')

    pathFilters = getPathFilters(channel);
    [offset, ~] = nrPerfectTimingEstimate(pathGains, pathFilters);
    hest = nrPerfectChannelEstimate(pathGains, pathFilters, NRB, SCS, nSlot, offset, sampleTimes);

    subplot(2, 2, 1);
    z1 = pow2db(abs(hest(:, :, 1, 1)).^2);
    surf(z1);
    Zmax1 = max(z1, [], 'all');
    shading('flat');
    xlabel('OFDM Symbols'); ylabel('Subcarriers'); zlabel('Magnitude Squared (dB)');
    title('Channel Magnitude Response (1^{st} tx - 1^{st} rx antenna)');
    disp("UE1 maximum Magnitude Response: " + Zmax1);

    [wbs, wue, ~] = getBeamformingWeights(hest, nLayers, scOffset, noRBs);

    % Plot UE radiation pattern
    ueSite.Antenna = clone(channel.ReceiveAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
    ueSite.Antenna.Taper = wue;
    pattern(ueSite, fc, "Size", 4);

    % Plot BS radiation pattern
    bsSite.Antenna = clone(channel.TransmitAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
    bsSite.Antenna.Taper = wbs;
    pattern(bsSite, fc, "Size", 5);
catch
    warning('UE1 cannot be ray traced');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ueSite2 = rxsite("Name", "UE2", ...
"Latitude", uePosition2(1), "Longitude", uePosition2(2), ...
    "AntennaHeight", 1, ... % in m
    "AntennaAngle", ueArrayOrientation2(1:2));

ueSite2.show();
rays2 = raytrace(bsSite, ueSite2, pm, "Type", "pathloss");

try
    plot(rays2{1})
    pathToAs2 = [rays2{1}.PropagationDelay] - min([rays2{1}.PropagationDelay]); % Time of arrival of each ray (normalized to 0 sec)
    avgPathGains2 =- [rays2{1}.PathLoss]; % Average path gains of each ray
    pathAoDs2 = [rays2{1}.AngleOfDeparture]; % AoD of each ray
    pathAoAs2 = [rays2{1}.AngleOfArrival]; % AoA of each ray
    isLOS2 = any([rays2{1}.LineOfSight]);

    channel2 = nrCDLChannel;
    channel2.DelayProfile = 'Custom';
    channel2.PathDelays = pathToAs2;
    channel2.AveragePathGains = avgPathGains2;
    channel2.AnglesAoD = pathAoDs2(1, :); % azimuth of departure
    channel2.AnglesZoD = 90 - pathAoDs2(2, :); % channel uses zenith angle, rays use elevation
    channel2.AnglesAoA = pathAoAs2(1, :); % azimuth of arrival
    channel2.AnglesZoA = 90 - pathAoAs2(2, :); % channel uses zenith angle, rays use elevation
    channel2.HasLOSCluster = isLOS2;
    channel2.CarrierFrequency = fc;
    channel2.NormalizeChannelOutputs = false; % do not normalize by the number of receive antennas, this would change the receive power
    channel2.NormalizePathGains = false; % set to false to retain the path gains

    channel2.SampleRate = waveformInfo.SampleRate;
    channel2.ReceiveAntennaArray = ueArray;
    channel2.ReceiveArrayOrientation = [ueArrayOrientation2(1); (-1) * ueArrayOrientation2(2); 0]; % the (-1) converts elevation to downtilt
    channel2.TransmitAntennaArray = bsArray;
    channel2.TransmitArrayOrientation = [bsArrayOrientation(1); (-1) * bsArrayOrientation(2); 0]; % the (-1) converts elevation to downtilt
    channel2.ChannelFiltering = false;

    [pathGains2, sampleTimes2] = channel2();
    pg2 = permute(pathGains2, [2 1 3 4]); % first dimension is the number of paths

    if isLOS2
        pg2 = [sum(pg2(1:2, :, :, :)); pg2(3:end, :, :, :)];
    end

    pg2 = abs(pg2).^2;

    pathFilters2 = getPathFilters(channel2);
    [offset2, ~] = nrPerfectTimingEstimate(pathGains2, pathFilters2);
    hest2 = nrPerfectChannelEstimate(pathGains2, pathFilters2, NRB, SCS, nSlot, offset2, sampleTimes2);
    [wbs2, wue2, ~] = getBeamformingWeights(hest2, nLayers, scOffset, noRBs);

    subplot(2, 2, 2);
    z2 = pow2db(abs(hest2(:, :, 1, 1)).^2);
    surf(z2);
    Zmax2 = max(z2, [], 'all');
    shading('flat');
    xlabel('OFDM Symbols'); ylabel('Subcarriers'); zlabel('Magnitude Squared (dB)');
    title('Channel Magnitude Response (1^{st} tx - 2^{nd} rx antenna)');
    disp("UE2 maximum Magnitude Response: " + Zmax2);

    ueSite2.Antenna = clone(channel2.ReceiveAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
    ueSite2.Antenna.Taper = wue2;
    pattern(ueSite2, fc, "Size", 4);

    %bsSite.Antenna = clone(channel2.TransmitAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
    bsSite.Antenna.Taper = wbs2;
    pattern(bsSite, fc, "Size", 5);
catch
    warning('UE2 cannot be ray traced');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ueSite3 = rxsite("Name", "UE3", ...
"Latitude", uePosition3(1), "Longitude", uePosition3(2), ...
    "AntennaHeight", 1, ... % in m
    "AntennaAngle", ueArrayOrientation3(1:2));

ueSite3.show();
rays3 = raytrace(bsSite, ueSite3, pm, "Type", "pathloss");

try
    plot(rays3{1})
    pathToAs3 = [rays3{1}.PropagationDelay] - min([rays3{1}.PropagationDelay]); % Time of arrival of each ray (normalized to 0 sec)
    avgPathGains3 =- [rays3{1}.PathLoss]; % Average path gains of each ray
    pathAoDs3 = [rays3{1}.AngleOfDeparture]; % AoD of each ray
    pathAoAs3 = [rays3{1}.AngleOfArrival]; % AoA of each ray
    isLOS3 = any([rays3{1}.LineOfSight]);

    channel3 = nrCDLChannel;
    channel3.DelayProfile = 'Custom';
    channel3.PathDelays = pathToAs3;
    channel3.AveragePathGains = avgPathGains3;
    channel3.AnglesAoD = pathAoDs3(1, :); % azimuth of departure
    channel3.AnglesZoD = 90 - pathAoDs3(2, :); % channel uses zenith angle, rays use elevation
    channel3.AnglesAoA = pathAoAs3(1, :); % azimuth of arrival
    channel3.AnglesZoA = 90 - pathAoAs3(2, :); % channel uses zenith angle, rays use elevation
    channel3.HasLOSCluster = isLOS3;
    channel3.CarrierFrequency = fc;
    channel3.NormalizeChannelOutputs = false; % do not normalize by the number of receive antennas, this would change the receive power
    channel3.NormalizePathGains = false; % set to false to retain the path gains

    channel3.SampleRate = waveformInfo.SampleRate;
    channel3.ReceiveAntennaArray = ueArray;
    channel3.ReceiveArrayOrientation = [ueArrayOrientation3(1); (-1) * ueArrayOrientation3(2); 0]; % the (-1) converts elevation to downtilt
    channel3.TransmitAntennaArray = bsArray;
    channel3.TransmitArrayOrientation = [bsArrayOrientation(1); (-1) * bsArrayOrientation(2); 0]; % the (-1) converts elevation to downtilt
    channel3.ChannelFiltering = false;

    [pathGains3, sampleTimes3] = channel3();
    pg3 = permute(pathGains3, [2 1 3 4]); % first dimension is the number of paths

    if isLOS3
        pg3 = [sum(pg3(1:2, :, :, :)); pg3(3:end, :, :, :)];
    end

    pg3 = abs(pg3).^2;

    pathFilters3 = getPathFilters(channel3);
    [offset3, ~] = nrPerfectTimingEstimate(pathGains3, pathFilters3);
    hest3 = nrPerfectChannelEstimate(pathGains3, pathFilters3, NRB, SCS, nSlot, offset3, sampleTimes3);
    [wbs3, wue3, ~] = getBeamformingWeights(hest3, nLayers, scOffset, noRBs);

    subplot(2, 2, 3);
    z3 = pow2db(abs(hest3(:, :, 1, 1)).^2);
    surf(z3);
    Zmax3 = max(z3, [], 'all');
    shading('flat');
    xlabel('OFDM Symbols'); ylabel('Subcarriers'); zlabel('Magnitude Squared (dB)');
    title('Channel Magnitude Response (1^{st} tx - 3^{rd} rx antenna)');
    disp("UE3 maximum Magnitude Response: " + Zmax3);
    ueSite3.Antenna = clone(channel3.ReceiveAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
    ueSite3.Antenna.Taper = wue3;
    pattern(ueSite3, fc, "Size", 4);

    bsSite.Antenna = clone(channel3.TransmitAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
    bsSite.Antenna.Taper = wbs3;
    pattern(bsSite, fc, "Size", 5);
catch
    warning('UE3 cannot be ray traced');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ueSite4 = rxsite("Name", "UE4", ...
"Latitude", uePosition4(1), "Longitude", uePosition4(2), ...
    "AntennaHeight", 1, ... % in m
    "AntennaAngle", ueArrayOrientation4(1:2));

ueSite4.show();
rays4 = raytrace(bsSite, ueSite4, pm, "Type", "pathloss");

try
    plot(rays4{1})
    pathToAs4 = [rays4{1}.PropagationDelay] - min([rays4{1}.PropagationDelay]); % Time of arrival of each ray (normalized to 0 sec)
    avgpathGains4 =- [rays4{1}.PathLoss]; % Average path gains of each ray
    pathAoDs4 = [rays4{1}.AngleOfDeparture]; % AoD of each ray
    pathAoAs4 = [rays4{1}.AngleOfArrival]; % AoA of each ray
    isLOS4 = any([rays4{1}.LineOfSight]);

    rays4 = nrCDLChannel;
    rays4.DelayProfile = 'Custom';
    rays4.PathDelays = pathToAs4;
    rays4.AveragePathGains = avgpathGains4;
    rays4.AnglesAoD = pathAoDs4(1, :); % azimuth of departure
    rays4.AnglesZoD = 90 - pathAoDs4(2, :); % channel uses zenith angle, rays use elevation
    rays4.AnglesAoA = pathAoAs4(1, :); % azimuth of arrival
    rays4.AnglesZoA = 90 - pathAoAs4(2, :); % channel uses zenith angle, rays use elevation
    rays4.HasLOSCluster = isLOS4;
    rays4.CarrierFrequency = fc;
    rays4.NormalizeChannelOutputs = false; % do not normalize by the number of receive antennas, this would change the receive power
    rays4.NormalizePathGains = false; % set to false to retain the path gains

    rays4.SampleRate = waveformInfo.SampleRate;
    rays4.ReceiveAntennaArray = ueArray;
    rays4.ReceiveArrayOrientation = [ueArrayOrientation4(1); (-1) * ueArrayOrientation4(2); 0]; % the (-1) converts elevation to downtilt
    rays4.TransmitAntennaArray = bsArray;
    rays4.TransmitArrayOrientation = [bsArrayOrientation(1); (-1) * bsArrayOrientation(2); 0]; % the (-1) converts elevation to downtilt
    rays4.ChannelFiltering = false;

    [pathGains4, sampleTimes4] = rays4();
    pg4 = permute(pathGains4, [2 1 3 4]); % first dimension is the number of paths

    if isLOS4
        pg4 = [sum(pg4(1:2, :, :, :)); pg4(3:end, :, :, :)];
    end

    pg4 = abs(pg4).^2;

    pathFilters4 = getPathFilters(rays4);
    [offset4, ~] = nrPerfectTimingEstimate(pathGains4, pathFilters4);
    hest4 = nrPerfectChannelEstimate(pathGains4, pathFilters4, NRB, SCS, nSlot, offset4, sampleTimes4);
    [wbs4, wue4, ~] = getBeamformingWeights(hest4, nLayers, scOffset, noRBs);

    subplot(2, 2, 4);
    z4 = pow2db(abs(hest4(:, :, 1, 1)).^2);
    surf(z4);
    Zmax4 = max(z4, [], 'all');
    shading('flat');
    xlabel('OFDM Symbols'); ylabel('Subcarriers'); zlabel('Magnitude Squared (dB)');
    title('Channel Magnitude Response (1^{st} tx - 4^{th} rx antenna)');
    disp("UE4 maximum Magnitude Response: " + Zmax4);
    ueSite4.Antenna = clone(rays4.ReceiveAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
    ueSite4.Antenna.Taper = wue4;
    pattern(ueSite4, fc, "Size", 4);

    bsSite.Antenna = clone(rays4.TransmitAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
    bsSite.Antenna.Taper = wbs4;
    pattern(bsSite, fc, "Size", 5);
catch
    warning('UE4 cannot be ray traced');
end

function [wtx, wrx, D] = getBeamformingWeights(hEst, nLayers, scOffset, noRBs)
    % Get beamforming weights given a channel matrix hEst and the number of
    % layers nLayers. One set of weights is provided for the whole bandwidth.
    % The beamforming weights are calculated using singular value (SVD)
    % decomposition.
    %
    % Only part of the channel estimate is used to get the weights, this is
    % indicated by an offset SCOFFSET (offset from the first subcarrier) and a
    % width in RBs (NORBS).

    % Average channel estimate
    [~, ~, R, P] = size(hEst);
    %H = permute(mean(reshape(hEst,[],R,P)),[2 3 1]);

    scNo = scOffset + 1;
    hEst = hEst(scNo:scNo + (12 * noRBs - 1), :, :, :);
    H = permute(mean(reshape(hEst, [], R, P)), [2 3 1]);

    % SVD decomposition
    [U, D, V] = svd(H);
    wtx = V(:, 1:nLayers).';
    wrx = U(:, 1:nLayers)';
end

%%
function [Carrier] = redcap_CarrierConfig(Carrier)
    %
    %
    %
    %
    %

    % The number of OFDM symbols in a slot depends on the cyclic prefix
    if strcmpi(Carrier.CyclicPrefix, 'normal')

        Carrier.SymbolsPerSlot = 14;
    else

        Carrier.SymbolsPerSlot = 12;
    end

    % The number of slots in a subframe depends on the subcarrier spacing
    Carrier.SlotsPerSubframe = double(Carrier.SubcarrierSpacing) / 15;

    % The number of slots in a frame depends on the subcarrier spacing
    Carrier.SlotsPerFrame = 10 * (double(Carrier.SubcarrierSpacing) / 15);

end

function info = redcap_OFDMInfo(varargin)
    %
    %
    %
    %   SubcarrierSpacing - Subcarrier spacing in kHz (15, 30, 60, 120, 240)
    %   CyclicPrefix      - Cyclic prefix ('normal', 'extended')
    %   NSizeGrid         - Number of resource blocks in carrier resource grid
    %                       (1...275)
    %
    %   INFO is a structure containing the fields:
    %
    %   Nfft                - Number of IFFT points used in the OFDM modulator
    %   SampleRate          - Sample rate of the OFDM modulated waveform
    %   CyclicPrefixLengths - Cyclic prefix length (in samples) of each OFDM
    %                         symbol in a subframe, starting at slot 0
    %   SymbolLengths       - Total length (in samples) of each OFDM symbol in
    %                         a subframe, including the cyclic prefix and
    %                         starting at slot 0
    %   Windowing           - Number of time-domain samples over which
    %                         windowing and overlapping of OFDM symbols is
    %                         applied
    %   SymbolPhases        - Phase precompensation applied for each OFDM
    %                         symbol due to the phase term per OFDM symbol in
    %                         TS 38.211 Section 5.4. <a
    %                         href="matlab:help('nrOFDMModulate')"
    %                         >nrOFDMModulate</a> applies
    %                         this precompensation during modulation and
    %                         <a href="matlab:help('nrOFDMDemodulate')"
    %                         >nrOFDMDemodulate</a> performs decompensation
    %                         during demodulation.
    %   SymbolsPerSlot      - Number of OFDM symbols in a slot
    %   SlotsPerSubframe    - Number of slots in a 1 ms subframe
    %   SlotsPerFrame       - Number of slots in a 10 ms frame
    %
    %   Note that the number of samples in the INFO.CyclicPrefixLengths,
    %   INFO.SymbolLengths, and INFO.Windowing fields apply to the sample rate
    %   of the IFFT of size INFO.Nfft used during OFDM symbol construction.
    %   This may be different from the sample rate of the waveform in the case
    %   that the 'SampleRate' NAME,VALUE pair below is specified. Note also
    %   that the IFFT size can be specified using the 'Nfft' NAME,VALUE pair.
    %
    %   INFO = nrOFDMInfo(NRB,SCS) provides dimensional information related to
    %   OFDM modulation for NRB resource blocks with subcarrier spacing SCS.
    %
    %   NRB is the number of resource blocks (1...275).
    %
    %   SCS is the subcarrier spacing in kHz (15, 30, 60, 120, 240).
    %
    %   INFO = nrOFDMInfo(...,NAME,VALUE) specifies additional options as
    %   NAME,VALUE pairs to allow control over the OFDM modulation:
    %
    %   CyclicPrefix        - Cyclic prefix ('normal' (default), 'extended').
    %                         This option is only applicable for function
    %                         syntaxes not using nrCarrierConfig
    %   Nfft                - Desired number of IFFT points to use in the OFDM
    %                         modulator. If absent or set to [], a default
    %                         value is selected based on other parameters, see
    %                         <a href="matlab: doc('nrOFDMModulate')"
    %                         >nrOFDMModulate</a> for details
    %   SampleRate          - Desired sample rate of the OFDM modulated
    %                         waveform. If absent or set to [], the default
    %                         value is SampleRate = Nfft * SCS. If required,
    %                         the OFDM modulated waveform is resampled to this
    %                         sample rate after OFDM symbol construction, using
    %                         an IFFT of size INFO.Nfft
    %   Windowing           - Number of time-domain samples over which
    %                         windowing and overlapping of OFDM symbols is
    %                         applied. If absent or set to [], a default value
    %                         is selected based on other parameters, see
    %                         <a href="matlab: doc('nrOFDMModulate')"
    %                         >nrOFDMModulate</a> for details
    %   CarrierFrequency    - Carrier frequency (in Hz) to calculate the phase
    %                         precompensation applied for each OFDM symbol
    %                         (denoted f_0 in TS 38.211 Section 5.4). Default
    %                         is 0
    %
    % version           author               data                 note
    % V1.2               sxj                     9/9         更改不以类输入时候的CP的参数赋�?

    narginchk(1, 12);

    % Validate inputs and get OFDM information
    internalinfo = validateInputs(varargin{:});

    % Create output structure
    info = nr5g.internal.OFDMInfoOutput(internalinfo);

end

% Validate inputs
function info = validateInputs(varargin)

    coder.extrinsic('nr5g.internal.parseOptions');

    fcnName = 'nrOFDMInfo';

    isCarrierSyntax = isa(varargin{1}, 'nrCarrierConfig');

    if (isCarrierSyntax) % nrOFDMInfo(CARRIER,...)

        % Validate carrier input type
        carrier = varargin{1};
        validateattributes(carrier, {'nrCarrierConfig'}, {'scalar'}, fcnName, 'Carrier specific configuration object');

        % Parse options
        optNames = {'Nfft', 'SampleRate', 'Windowing', 'CarrierFrequency'};
        opts = coder.const(nr5g.internal.parseOptions([fcnName '(carrier,...'], optNames, varargin{2:end}));

        % Get OFDM information
        info = nr5g.internal.OFDMInfo(carrier, opts);

    else % nrOFDMInfo(NRB,SCS,...)

        % Validate NRB
        NRB = varargin{1};
        validateattributes(NRB, {'numeric'}, {'real', 'integer', 'scalar', '>=', 1, '<=', 275}, fcnName, 'NRB');

        % Validate subcarrier spacing
        SCS = varargin{2};
        validateattributes(SCS, {'numeric'}, {'real', 'integer', 'scalar'}, fcnName, 'SCS');
        validSCS = [15 30 60 120 240];
        coder.internal.errorIf(~any(SCS == validSCS), 'nr5g:nrOFDMInfo:InvalidSCS', SCS, num2str(validSCS));

        % Parse options and get cyclic prefix length
        %         optNames = {'CyclicPrefix','Nfft','SampleRate','Windowing','CarrierFrequency'};
        %         opts = coder.const(nr5g.internal.parseOptions(fcnName,optNames,varargin{3:end}));

        optNames = {'Nfft', 'SampleRate', 'Windowing', 'CarrierFrequency'};
        opts = coder.const(nr5g.internal.parseOptions(fcnName, optNames, varargin{4:end}));

        ECP = strcmpi(varargin{3}, 'extended');

        % Get OFDM information
        %info = coder.const(feval('nr5g.internal.OFDMInfo',NRB,SCS,ECP,opts));
        info = nr5g.internal.OFDMInfo(NRB, SCS, ECP, opts);

    end

end
