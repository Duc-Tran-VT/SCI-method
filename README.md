# SCI-method

All data and code in this repository are associated with the following paper:

Tran, D., Jacquet, M., Pearson, S., van Prooijen, B., and Verney, R. (2024). Estimation of
mud and sand fractions and total concentration from coupled optical-acoustic sensors.
(https://zenodo.org/doi/10.5281/zenodo.12819148).

Laboratory data used in this study obtained from a series of experiments conducted in the DEXMES facility in IFREMER.

* Note 1. The data is raw, uncalibrated signal. More detailed information is in Section 2.2 Data processing.

* Note 2. Files' names and what does it do

All.mat: Data from Q_set. 

S1S2.mat: Data from V_set

ABS.m: Processing data related to the sensor ABS and all optical sensors

ADV.m: Processing data related to the sensor ADV and all optical sensors

AQUA.m: Processing data related to the sensor AQUAscat and all optical sensors

* Note 3. Notation used in the data file and in text of the manuscript:

abs: LISST_ABS 8 Mhz (Laser In-Situ Scattering and Transmissiometery - Acoustic Backscatter Sensor) (in manuscript read: A_8, unit: mg/L)

adv: Nortek Vector Acoustic Doppler Velocimeter (in manuscript read: A_6, unit: dB)

aqua4: AQUAscat-1000R, transducer 4 Mhz (in manuscript read: A_4, unit: count)

aqua3: AQUAscat-1000R, transducer 2 Mhz (in manuscript read: A_2, unit: count)

aqua2: AQUAscat-1000R, transducer 1 Mhz (in manuscript read: A_1, unit: count)

aqua1: AQUAscat-1000R, transducer 0.5 Mhz (in manuscript read: A_0.5, unit: count)

hydro4: HydroScat-4, channel 852 nm (in manuscript read: O_852, unit: m-1)

hydro3: HydroScat-4, channel 620 nm (in manuscript read: O_620, unit: m-1)

hydro2: HydroScat-4, channel 532 nm (in manuscript read: O_532, unit: m-1)

hydro1: HydroScat-4, channel 420 nm (in manuscript read: O_420, unit: m-1)

wet: Wetlabs_FLNTU 700 nm (in manuscript read: O_700, unit: NTU)

C: concentration obtained from filtered water sample (in manuscript read: C, unit: mg/L)
