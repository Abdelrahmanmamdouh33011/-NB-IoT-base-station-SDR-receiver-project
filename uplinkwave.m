function waveform = uplinkwave()

%% NB-IoT Uplink Waveform Generation
% This example shows how to generate LTE-Advanced Pro Release 13 Narrowband
% IoT (NB-IoT) uplink waveforms consisting of the Narrowband Physical
% Uplink Shared Channel (NPUSCH) and the associated demodulation reference
% signals for test and measurement applications using the LTE Toolbox(TM).

% Copyright 2018 The MathWorks, Inc.

%% Introduction 
% 3GPP introduced a new air interface, Narrowband IoT (NB-IoT) optimized
% for low data rate machine type communications in LTE-Advanced Pro Release
% 13. NB-IoT provides cost and power efficiency improvements as it avoids
% the need for complex signaling overhead required for LTE based systems.
%
% The LTE Toolbox can be used to generate standard compliant NB-IoT uplink
% complex baseband waveforms representing the 180kHz narrowband carrier
% suitable for test and measurement applications. The LTE Toolbox supports
% all the NB-IoT modes of operation described below - standalone, guardband
% and in-band.
%
% * Standalone: NB-IoT carrier deployed outside the LTE spectrum, e.g.
% the spectrum used for GSM or satellite communications
% * Guardband: NB-IoT carrier deployed in the guardband between two LTE
% carriers
% * In-band: NB-IoT carrier deployed in resource blocks of an LTE carrier
%
% The NB-IoT uplink consists of the following physical layer channels and
% signals:
% 
% * Narrowband demodulation reference signal (DM-RS)
% * Narrowband physical uplink shared channel (NPUSCH)
% * Narrowband physical random access channel (NPRACH)
%
% In this example, we demonstrate the NB-IoT uplink resource element (RE)
% grid and waveform generation consisting of the NPUSCH and DM-RS signals.
% The sections below introduce these physical signals and channels that
% form the grid along with key concepts including subframe repetition,
% logical and transport channel mappings, and the corresponding grids for
% the different configurations.
%
% The example outputs the complex baseband waveform along with the
% populated grid containing the NPUSCH and DM-RS signals. The waveform can
% be used for a range of applications from RF testing to simulation of
% receiver implementations.

%% NPUSCH Allocation
% This section provides an overall description of how the NPUSCH gets
% mapped into the NB-IoT uplink slots. 
% 
% The NPUSCH can carry the uplink shared channel (UL-SCH) or the uplink
% control information according to the two formats: 
% 
% * NPUSCH format 1, used to carry the uplink shared channel (UL-SCH) 
% * NPUSCH format 2, used to carry uplink control information
% 
% The NPUSCH is transmitted on one or more resource units and each of these
% resource units are repeated up to 128 times to improve transmission
% reliability and coverage without compromising on the low power and low
% complexity requirements to meet ultra low end IoT use cases.
%
% The smallest mapping unit for the NPUSCH is a resource unit. It is
% defined as 7 * |NslotsUL| consecutive SC-FDMA symbols in the time domain
% and |NscRU| consecutive subcarriers in the frequency domain, where
% |NslotsUL| and |NscRU| are defined in TS 36.211 Table 10.1.2.3-1 [ <#9 1>
% ]. The NB-IoT UL-SCH may carry Common Control Channel (CCCH), Dedicated
% Control Channel (DCCH) or Dedicated Traffic Channel (DTCH) and maps on to
% the NPUSCH physical channel (TS 36.300 section 6.1.3.1 and section 5.3.1a
% [ <#9 6> ]). The NPUSCH can be mapped to one or more than one resource
% units, |NRU| as defined by TS 36.211 Section 10.1.3.6 [ <#9 1> ] and each
% resource unit can be transmitted |Nrep| times.
% 
% The examples in the figure shows the repetition pattern with |NRep| = 4.
% The total duration for transmitting a block of data is given by |NRU| *
% |NULSlots| * |MidenticalNPUSCH| as specified in TS 36.211 Section
% 10.1.3.6 [ <#9 1> ]. For the first case shown below, each transport block
% is transmitted over |NRU| = 2 and each of these |NRU| contains two UL
% slots indicated by |NULSlots|. After mapping to |Nslots|, these slots
% will be repeated |MidenticalNPUSCH| = 2 (assuming NscRU > 1) times. In
% the second case, we assume that |NscRU| is 1 and hence |MidenticalNPUSCH|
% = 1. This, combined with |Nslots| = 1 results in the transmission pattern
% where each block is transmitted without internal repetitions. In all
% cases, the scrambling sequence is reset at the start of the codeword
% transmission or retransmission (see TS 36.211 Section 10.1.3.1 [ <#9 1>
% ]). The detailed specification of the repetition scheme can be found in
% TS 36.211 10.1.3 [ <#9 1> ].
%
% <<nbwaveformgen_npusch.png>>

%%  NB-IoT Uplink Slot Grid
% In addition to the slot allocation described above, this section further
% explains the RE allocation in a slot. The grid consists of one or more
% frames containing NPUSCH and corresponding DM-RS.
%
% * _DM-RS_: The DM-RS is transmitted in every NPUSCH slot with the same
% bandwidth as the associated NPUSCH. The reference signals depend on the
% number of subcarriers |NscRU|, the narrowband cell ID |NNcellID| and the
% NPUSCH format |NPUSCHFormat|. The RE positions depend on the NPUSCH
% format and subcarrier spacing. For NPUSCH format 1 with subcarrier
% spacing of 3.75kHz, the DM-RS is transmitted on symbol 4 and with
% subcarrier spacing of 15kHz, the DM-RS is transmitted on symbol 3. For
% NPUSCH format 2 with subcarrier spacing of 3.75kHz, the DM-RS is
% transmitted on symbols 0,1,2 and with subcarrier spacing of 15kHz, the
% DM-RS is transmitted on symbols 2,3 and 4 in the slot.
% * _NPUSCH_: The NPUSCH supports single tone bandwidth in addition to the
% multitone (12 subcarriers) bandwidth. Single tone transmissions can use
% either the 15kHz or 3.75kHz subcarrier spacing whereas the multitone
% transmissions use the 15kHz subcarrier spacing. This means that the slot
% duration for 15kHz mode is 0.5ms and for 3.75kHz the slot duration is
% 2ms. The scrambling sequence is initialized in the first slot of the
% transmission of the codeword. If there are repetitions enabled, then the
% scrambling sequence is reinitialized after every |MidenticalNPUSCH|
% transmission of the codeword as described in TS 36.211 Section 10.1.3.1[
% <#9 1> ]. The codewords are BPSK/QPSK modulated on to a single layer and
% precoded before mapping to one or more resource units. All the resource
% elements other than those used for demodulation reference signals are
% used for NPUSCH transmission. If higher layer signaling
% (|npusch-AllSymbols| as described in TS 36.211 Section 10.1.3.6 [ <#9 1>
% ]) indicates the presence of the SRS symbol, these symbols are counted in
% the NPUSCH mapping, but not used for the transmission of the NPUSCH (i.e.
% these NPUSCH positions are punctured by SRS).

%% NPUSCH Configuration
% In this section we configure the parameters required for NPUSCH
% generation. The UE uses the combination of MCS (modulation and coding
% scheme) and resource assignment signaled via the DCI to determine the
% transport block size from the set defined in TS 36.213 Table 16.5.1.2-2 [
% <#9 3> ] to use for NPUSCH transmission. In this example, this is
% specified via the parameter |tbs| and the duration of the generated
% waveform is controlled via the |totNumBlks| parameter.

tbs = 144;                        % The transport block size
totNumBlks = 1;                   % Number of simulated transport blocks

ue = struct();                    % Initialize the UE structure
ue.NBSubcarrierSpacing = '15kHz'; % 3.75kHz, 15kHz
ue.NNCellID = 0;                  % Narrowband cell identity

chs = struct();
% NPUSCH carries data or control information
chs.NPUSCHFormat = 'Data'; % Payload type (Data or Control)
% The number of subcarriers used for NPUSCH 'NscRU' depends on the NPUSCH
% format and subcarrier spacing 'NBSubcarrierSpacing' as shown in TS 36.211
% Table 10.1.2.3-1. There are 1,3,6 or 12 continuous subcarriers for NPUSCH
chs.NBULSubcarrierSet = 0:11;        % Range is 0-11 (15kHz); 0-47 (3.75kHz)
chs.NRUsc = length(chs.NBULSubcarrierSet);
chs.CyclicShift = 0;   % Cyclic shift required when NRUsc = 3 or 6
chs.RNTI = 0;          % RNTI value
chs.NLayers = 1;       % Number of layers 
chs.NRU = 2;           % Number of resource units
chs.NRep = 4;          % Number of repetitions of the NPUSCH
% The symbol modulation depends on the NPUSCH format and NscRU as
% given by TS 36.211 Table 10.1.3.2-1
chs.Modulation = 'QPSK';
rvDCI = 0;             % RV offset signaled via DCI (See 36.213 16.5.1.2)
if strcmpi(chs.NPUSCHFormat,'Control')
    if length(chs.NBULSubcarrierSet) ~= 1
        error('The number of useful subcarriers must be 1 for control information');
    end
end

% Specify the NPUSCH and DM-RS power scaling in dB for plot visualization
chs.NPUSCHPower = 30; 
chs.NPUSCHDRSPower = 34;

%%
% For DM-RS signals in NPUSCH format 1, sequence-group hopping can be
% enabled or disabled by the higher layer cell-specific parameter
% |groupHoppingEnabled|. Sequence-group hopping for a particular UE can be
% disabled through the higher layer parameter |groupHoppingDisabled| as
% described in TS 36.211 Section 10.1.4.1.3 [ <#9 1> ]. In this example, we
% use the |SeqGroupHopping| parameter to enable or disable sequence-group
% hopping
chs.SeqGroupHopping = true; % Enable/Disable Sequence-Group Hopping for UE
chs.Deltass = 0;            % Higher-layer parameter groupAssignmentNPUSCH

% Get number of time slots in a resource unit NULSlots
if strcmpi(chs.NPUSCHFormat,'Data')
    if strcmpi(ue.NBSubcarrierSpacing,'3.75kHz')
        NULSlots = 16;
    elseif strcmpi(ue.NBSubcarrierSpacing,'15kHz')
        if chs.NRUsc == 12
            NULSlots = 2;
        elseif chs.NRUsc == 6
            NULSlots = 4;
        elseif chs.NRUsc == 3
            NULSlots = 8;
        elseif chs.NRUsc == 1
            NULSlots = 16;
        else
            error('Invalid number of subcarriers, should be one of 1,3,6,12');
        end
    else
        error('Invalid subcarrier spacing, should be either 3.75kHz or 15kHz');
    end
elseif strcmpi(chs.NPUSCHFormat,'Control')
    NULSlots = 4;
else
    error('Invalid NPUSCH Format (%s), should be either Data or Control',chs.NPUSCHFormat);
end
chs.NULSlots = NULSlots;
NSlotsPerBundle = chs.NRU*chs.NULSlots*chs.NRep; % Number of slots in a codeword bundle
TotNSlots = totNumBlks*NSlotsPerBundle;   % Total number of simulated slots

%% NB-IoT Uplink Waveform Generation
% In this section, we create the resource grid populated with the NPUSCH
% and the corresponding demodulation reference signals. This grid is then
% SC-FDMA modulated to generate the time domain waveform.

% Initialize the random generator to default state 
rng('default');

% Get the slot grid and number of slots per frame
if strcmpi(ue.NBSubcarrierSpacing,'15kHz')
    slotGrid = zeros(12, 7);
    NSlotsPerFrame = 20; % Slots 0...19 
else
    slotGrid = zeros(48, 7);
    NSlotsPerFrame = 5;  % Slots 0...4
end

estate = [];    % NPUSCH encoder state, auto re-initialization in the function
trblk = [];     % Initialize the transport block
txgrid = [];    % Full grid initialization
statedmrs = []; % DM-RS state

% Display the number of slots being generated
fprintf('\nGenerating %d slots corresponding to %d transport block(s)\n',TotNSlots,totNumBlks);
for slotIdx = 0+(0:TotNSlots-1)
    % Calculate the frame number and slot number within the frame
    ue.NFrame = fix(slotIdx/NSlotsPerFrame);
    ue.NSlot = mod(slotIdx,NSlotsPerFrame);

    if isempty(trblk)
       if strcmpi(chs.NPUSCHFormat,'Data')
           % UL-SCH encoding is done for the two RV values used for
           % transmitting the codewords. The RV sequence used is determined
           % from the rvDCI value signaled in the DCI and alternates
           % between 0 and 2 as given in TS 36.213 Section 16.5.1.2
                     
           % Define the transport block which will be encoded to create the
           % codewords for different RV
           trblk = randi([0 1],tbs,1);
            
           % Determine the coded transport block size
           [~, info] = hNPUSCHIndices(ue,chs);
           outblklen = info.G;
           % Create the codewords corresponding to the two RV values used
           % in the first and second block, this will be repeated till all
           % blocks are transmitted
           chs.RV = 2*mod(rvDCI+0,2); % RV for the first block
           cw = hNULSCH(chs,outblklen,trblk); % CRC and Turbo coding is repeated
           chs.RV = 2*mod(rvDCI+1,2); % RV for the second block
           cw = [cw hNULSCH(chs,outblklen,trblk)]; %#ok<AGROW> % CRC and Turbo coding is repeated
       else
           trblk = randi([0 1],1); % 1 bit ACK
           % For ACK, the same codeword is transmitted every block as
           % defined in TS 36.212 Section 6.3.3
           cw = hNULSCH(trblk);
       end
       blockIdx = 0; % First block to be transmitted
    end

    % NPUSCH encoding and mapping onto the slot grid
    [txsym, estate] = hNPUSCH(ue,chs,cw(:,mod(blockIdx,size(cw,2))+1),estate);
    % If a full block is transmitted, increment the clock counter so that
    % the correct codeword can be selected
    if estate.EndOfBlk
        blockIdx = blockIdx + 1;
    end
    % Map NPUSCH symbols in the grid of a slot
    indicesNPUSCH = hNPUSCHIndices(ue,chs);
    slotGrid(indicesNPUSCH) = txsym*db2mag(chs.NPUSCHPower);
    
    % Create DM-RS sequence and map to the slot grid
    [dmrs,statedmrs] = hNPUSCHDRS(ue,chs,statedmrs);
    indicesDMRS = hNPUSCHDRSIndices(ue,chs);
    slotGrid(indicesDMRS) = dmrs*db2mag(chs.NPUSCHDRSPower);

    % Concatenate this slot to the slot grid
    txgrid = [txgrid slotGrid]; %#ok<AGROW>
    
    % trblk err count and re-initialization
    if estate.EndOfTx
       % Re-initialize to enable the transmission of a new transport block
       trblk = [];
    end
    
end

% Perform SC-FDMA modulation to create time domain waveform
ue.CyclicPrefixUL = 'Normal'; % Normal cyclic prefix length for NB-IoT
[waveform,scfdmaInfo] = lteSCFDMAModulate(ue,chs,txgrid);

% %% Plot Transmitted Grid
% % Plot the populated grid and observe the NPUSCH and corresponding DM-RS.
% % The positions of the NPUSCH and DM-RS depends on the number of
% % subcarriers |chs.NRUsc| and the subcarriers used as specified by
% % |chs.NBULSubcarrierSet|. Note that the resource grid plot uses the power
% % levels of the PUSCH and the DM-RS to assign colors to the resource
% % elements.
% 
% % Create an image of overall resource grid
% figure
% image(abs(txgrid(:,:,1))); 
% axis xy;
% title(sprintf('NB-IoT Uplink RE Grid (NRep = %d, NRUsc = %d, NRU = %d)',chs.NRep,chs.NRUsc,chs.NRU)) 
% xlabel('OFDM symbols')
% ylabel('Subcarriers')
% % Create the legend box to indicate the channel/signal types associated with the REs
% reNames = {'NPUSCH';'DM-RS'};
% clevels = round(db2mag([chs.NPUSCHPower chs.NPUSCHDRSPower]));
% N = numel(reNames);
% L = line(ones(N),ones(N), 'LineWidth',8); % Generate lines                   
% % Index the color map
% cmap = colormap(gcf);
% % Set the colors according to cmap
% set(L,{'color'},mat2cell(cmap( min(1+clevels,length(cmap) ),:),ones(1,N),3));
% legend(reNames{:});

%% Appendix
% This example uses the helper functions:
%
% * <matlab:edit('hNULSCH.m') hNULSCH.m>
% * <matlab:edit('hNPUSCH.m') hNPUSCH.m>
% * <matlab:edit('hNPUSCHIndices.m') hNPUSCHIndices.m>
% * <matlab:edit('hNPUSCHDRS.m') hNPUSCHDRS.m>
% * <matlab:edit('hNPUSCHDRSIndices.m') hNPUSCHDRSIndices.m>

%% Selected Bibliography
% # 3GPP TS 36.211 "Physical channels and modulation"
% # 3GPP TS 36.212 "Multiplexing and channel coding"
% # 3GPP TS 36.213 "Physical layer procedures"
% # 3GPP TS 36.321 "Medium Access Control (MAC); Protocol specification"
% # 3GPP TS 36.331 "Radio Resource Control (RRC); Protocol specification"
% # 3GPP TS 36.300 "Overall description; Stage 2"
% # O. Liberg, M. Sundberg, Y.-P. Wang, J. Bergman and J. Sachs, Cellular Internet of Things: Technologies, Standards and Performance, Elsevier, 2018.

displayEndOfDemoMessage(mfilename)
end

