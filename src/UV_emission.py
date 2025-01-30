from pylab import *
import pandas as pd
import netCDF4 as nc
from scipy.interpolate import interp1d

#Absorption and Rayleigh scattering cross sections from VULCAN chemical kinetics program (Tsai et al., 2017)
H2_absdat = '/data4/ananyo/models/C3M/data/VULCAN/H2/H2_cross.csv'
He_absdat = '/data4/ananyo/models/C3M/data/VULCAN/He/He_cross.csv'
H2_scatdat = '/data4/ananyo/models/C3M/data/VULCAN/rayleigh/H2_scat.txt'
He_scatdat = '/data4/ananyo/models/C3M/data/VULCAN/rayleigh/He_scat.txt'
CH4_absdat = '/data4/ananyo/models/C3M/data/VULCAN/CH4/CH4_cross.csv'
C2H2_absdat = '/data4/ananyo/models/C3M/data/VULCAN/C2H2/C2H2_cross.csv'

data1 = pd.read_csv(H2_absdat, skiprows = 2, usecols = [0,1,2,3], names = ['wav', 'abs','diss','ion'])

data2 = pd.read_csv(He_absdat, skiprows = 2, usecols = [0,1,2,3], names = ['wav','abs','diss','ion'])

data3 = genfromtxt(H2_scatdat, skip_header=1)

data4 = genfromtxt(He_scatdat, skip_header=1)

data5 = pd.read_csv(CH4_absdat, skiprows = 2, usecols = [0,1,2,3], names = ['wav','abs','diss','ion'])

data6 = pd.read_csv(C2H2_absdat, skiprows = 2, usecols = [0,1,2,3], names = ['wav','abs','diss','ion'])

w1 = data1['wav'][:]*1e-9
ab1 = data1['abs'][:]*1e-4


w2 = data2['wav'][:]*1e-9
ab2 = data2['abs'][:]*1e-4

w5 = data5['wav'][:]*1e-9
ab5 = data5['abs'][:]*1e-4

w6 = data6['wav'][:]*1e-9
ab6 = data6['abs'][:]*1e-4

w3 = data3[:, 0]*1e-9
ab3 = data3[:,1]*1e-4

w4 = data4[:, 0]*1e-9
ab4 = data4[:,1]*1e-4

sol = genfromtxt('/home/ananyo/models/C3M/data/stellar/sun_spec.inp', skip_header=1)
w_sol = sol[:, 0]*1e-10
w_sol = linspace(50, 250, 100000)*1e-9

#Interpolate to solar grid
c1 = interp1d(w1, (ab1), fill_value=(0, 0), bounds_error=False)
c2 = interp1d(w2, (ab2), fill_value=(0, 0), bounds_error=False)
c3 = interp1d(w3, (ab3), fill_value=(0, 0), bounds_error=False)
c4 = interp1d(w4, (ab4), fill_value=(0, 0), bounds_error=False)
c5 = interp1d(w5, (ab5), fill_value=(0, 0), bounds_error=False)
c6 = interp1d(w6, (ab6), fill_value=(0, 0), bounds_error=False)

a1 = c1(w_sol)
a2 = c2(w_sol)
a3 = c3(w_sol)
a4 = c4(w_sol)
a5 = c5(w_sol)
a6 = c6(w_sol)

#Model atmosphere from Moses and Poppe 2017
#KINETICS model output
Pgrid = [6.708E+03, 5.976E+03, 5.305E+03, 4.689E+03, 4.122E+03, 3.706E+03, 3.414E+03, 3.137E+03,
 2.876E+03, 2.630E+03, 2.402E+03, 2.185E+03, 1.983E+03, 1.795E+03, 1.675E+03, 1.562E+03,
 1.455E+03, 1.352E+03, 1.255E+03, 1.162E+03, 1.074E+03, 9.901E+02, 9.113E+02, 8.370E+02,
 7.668E+02, 7.004E+02, 6.384E+02, 5.805E+02, 5.270E+02, 4.782E+02, 4.323E+02, 3.906E+02,
 3.515E+02, 3.160E+02, 2.831E+02, 2.531E+02, 2.251E+02, 1.999E+02, 1.770E+02, 1.563E+02,
 1.380E+02, 1.220E+02, 1.079E+02, 9.567E+01, 7.594E+01, 5.772E+01, 4.223E+01, 3.130E+01,
 2.340E+01, 1.763E+01, 1.275E+01, 9.287E+00, 6.818E+00, 5.039E+00, 3.751E+00, 2.803E+00,
 2.018E+00, 1.456E+00, 1.053E+00, 7.609E-01, 5.507E-01, 3.984E-01, 2.888E-01, 2.091E-01,
 1.517E-01, 1.100E-01, 7.977E-02, 5.787E-02, 4.206E-02, 3.056E-02, 2.222E-02, 1.619E-02,
 1.179E-02, 8.605E-03, 6.291E-03, 4.609E-03, 3.385E-03, 2.497E-03, 1.849E-03, 1.378E-03,
 9.982E-04, 7.319E-04, 5.446E-04, 4.006E-04, 2.945E-04, 2.200E-04, 1.653E-04, 1.228E-04,
 9.084E-05, 6.700E-05, 4.939E-05, 3.601E-05, 2.623E-05, 1.907E-05, 1.390E-05, 1.009E-05,
 7.311E-06, 5.285E-06, 3.806E-06, 2.743E-06, 1.978E-06, 1.425E-06, 1.025E-06, 7.390E-07,
 5.313E-07, 3.832E-07, 2.753E-07, 1.981E-07, 1.428E-07, 1.025E-07, 7.373E-08]

Hgrid = [-5.920E+01, -5.454E+01, -4.991E+01, -4.528E+01, -4.062E+01, -3.691E+01, -3.413E+01, -3.135E+01,
-2.856E+01, -2.576E+01, -2.299E+01, -2.020E+01, -1.741E+01, -1.464E+01, -1.278E+01, -1.091E+01,
-9.065E+00, -7.205E+00, -5.361E+00, -3.500E+00, -1.643E+00, 2.190E-01, 2.076E+00, 3.926E+00,
 5.784E+00, 7.645E+00, 9.498E+00, 1.136E+01, 1.322E+01, 1.507E+01, 1.694E+01, 1.879E+01,
 2.065E+01, 2.252E+01, 2.437E+01, 2.622E+01, 2.809E+01, 2.995E+01, 3.180E+01, 3.367E+01,
 3.552E+01, 3.737E+01, 3.924E+01, 4.110E+01, 4.480E+01, 4.945E+01, 5.502E+01, 6.059E+01,
 6.616E+01, 7.173E+01, 7.824E+01, 8.475E+01, 9.125E+01, 9.775E+01, 1.043E+02, 1.108E+02,
 1.182E+02, 1.257E+02, 1.331E+02, 1.406E+02, 1.481E+02, 1.555E+02, 1.630E+02, 1.705E+02,
 1.780E+02, 1.855E+02, 1.931E+02, 2.007E+02, 2.083E+02, 2.160E+02, 2.236E+02, 2.313E+02,
 2.391E+02, 2.469E+02, 2.547E+02, 2.626E+02, 2.706E+02, 2.786E+02, 2.867E+02, 2.948E+02,
 3.040E+02, 3.132E+02, 3.225E+02, 3.328E+02, 3.442E+02, 3.566E+02, 3.709E+02, 3.883E+02,
 4.086E+02, 4.317E+02, 4.577E+02, 4.875E+02, 5.202E+02, 5.556E+02, 5.930E+02, 6.331E+02,
 6.752E+02, 7.192E+02, 7.650E+02, 8.119E+02, 8.596E+02, 9.084E+02, 9.580E+02, 1.008E+03,
 1.058E+03, 1.109E+03, 1.161E+03, 1.212E+03, 1.264E+03, 1.316E+03, 1.369E+03]

Tgrid = [3.000E+02, 2.898E+02, 2.795E+02, 2.692E+02, 2.590E+02, 2.507E+02, 2.444E+02, 2.382E+02,
 2.318E+02, 2.255E+02, 2.194E+02, 2.129E+02, 2.064E+02, 2.002E+02, 1.960E+02, 1.917E+02,
 1.873E+02, 1.829E+02, 1.788E+02, 1.743E+02, 1.701E+02, 1.654E+02, 1.613E+02, 1.569E+02,
 1.524E+02, 1.476E+02, 1.443E+02, 1.417E+02, 1.390E+02, 1.366E+02, 1.340E+02, 1.313E+02,
 1.283E+02, 1.251E+02, 1.217E+02, 1.183E+02, 1.150E+02, 1.121E+02, 1.099E+02, 1.088E+02,
 1.087E+02, 1.100E+02, 1.119E+02, 1.143E+02, 1.201E+02, 1.269E+02, 1.328E+02, 1.375E+02,
 1.414E+02, 1.446E+02, 1.478E+02, 1.510E+02, 1.544E+02, 1.579E+02, 1.609E+02, 1.630E+02,
 1.644E+02, 1.650E+02, 1.653E+02, 1.654E+02, 1.654E+02, 1.654E+02, 1.654E+02, 1.654E+02,
 1.654E+02, 1.654E+02, 1.654E+02, 1.654E+02, 1.656E+02, 1.657E+02, 1.660E+02, 1.664E+02,
 1.668E+02, 1.675E+02, 1.683E+02, 1.694E+02, 1.709E+02, 1.728E+02, 1.754E+02, 1.789E+02,
 1.842E+02, 1.915E+02, 2.014E+02, 2.169E+02, 2.410E+02, 2.828E+02, 3.305E+02, 3.801E+02,
 4.327E+02, 4.865E+02, 5.395E+02, 5.925E+02, 6.419E+02, 6.870E+02, 7.267E+02, 7.618E+02,
 7.916E+02, 8.168E+02, 8.378E+02, 8.548E+02, 8.685E+02, 8.795E+02, 8.882E+02, 8.951E+02,
 9.006E+02, 9.048E+02, 9.082E+02, 9.107E+02, 9.128E+02, 9.144E+02, 9.157E+02]



H = [4.029E+02, 1.670E+03, 5.329E+03, 1.264E+04, 2.331E+04, 3.589E+04, 4.409E+04, 6.956E+04,
 1.052E+05, 1.216E+05, 1.504E+05, 1.745E+05, 2.078E+05, 2.423E+05, 2.714E+05, 3.019E+05,
 3.356E+05, 3.738E+05, 4.163E+05, 4.633E+05, 5.140E+05, 5.642E+05, 6.155E+05, 6.639E+05,
 7.117E+05, 7.607E+05, 8.324E+05, 9.253E+05, 1.035E+06, 1.164E+06, 1.316E+06, 1.486E+06,
 1.679E+06, 1.913E+06, 2.278E+06, 4.513E+06, 6.421E+06, 8.171E+06, 9.821E+06, 1.137E+07,
 1.274E+07, 1.382E+07, 1.456E+07, 1.484E+07, 1.393E+07, 1.164E+07, 9.512E+06, 7.972E+06,
 6.839E+06, 6.008E+06, 5.279E+06, 4.703E+06, 4.245E+06, 3.884E+06, 3.600E+06, 3.368E+06,
 3.148E+06, 2.973E+06, 2.846E+06, 2.771E+06, 2.749E+06, 2.785E+06, 2.895E+06, 3.111E+06,
 3.471E+06, 4.030E+06, 4.867E+06, 6.070E+06, 7.743E+06, 1.042E+07, 1.693E+07, 4.468E+07,
 1.563E+08, 4.481E+08, 9.344E+08, 1.504E+09, 2.007E+09, 2.335E+09, 2.462E+09, 2.401E+09,
 2.186E+09, 1.894E+09, 1.593E+09, 1.285E+09, 1.004E+09, 7.553E+08, 5.794E+08, 4.532E+08,
 3.593E+08, 2.882E+08, 2.340E+08, 1.904E+08, 1.562E+08, 1.290E+08, 1.073E+08, 8.940E+07,
 7.481E+07, 6.272E+07, 5.262E+07, 4.427E+07, 3.731E+07, 3.147E+07, 2.654E+07, 2.243E+07,
 1.893E+07, 1.602E+07, 1.352E+07, 1.143E+07, 9.665E+06, 8.158E+06, 6.894E+06]

He = [2.203E+19, 2.032E+19, 1.870E+19, 1.716E+19, 1.568E+19, 1.457E+19, 1.376E+19, 1.298E+19,
 1.223E+19, 1.149E+19, 1.078E+19, 1.011E+19, 9.464E+18, 8.832E+18, 8.422E+18, 8.027E+18,
 7.653E+18, 7.281E+18, 6.914E+18, 6.567E+18, 6.219E+18, 5.898E+18, 5.566E+18, 5.256E+18,
 4.957E+18, 4.676E+18, 4.359E+18, 4.036E+18, 3.736E+18, 3.449E+18, 3.178E+18, 2.931E+18,
 2.699E+18, 2.489E+18, 2.291E+18, 2.108E+18, 1.928E+18, 1.756E+18, 1.586E+18, 1.415E+18,
 1.249E+18, 1.091E+18, 9.485E+17, 8.234E+17, 6.218E+17, 4.469E+17, 3.122E+17, 2.232E+17,
 1.621E+17, 1.192E+17, 8.413E+16, 5.985E+16, 4.284E+16, 3.085E+16, 2.245E+16, 1.648E+16,
 1.169E+16, 8.338E+15, 5.969E+15, 4.270E+15, 3.058E+15, 2.188E+15, 1.567E+15, 1.120E+15,
 8.005E+14, 5.717E+14, 4.079E+14, 2.907E+14, 2.070E+14, 1.472E+14, 1.045E+14, 7.415E+13,
 5.249E+13, 3.711E+13, 2.621E+13, 1.848E+13, 1.301E+13, 9.154E+12, 6.426E+12, 4.504E+12,
 3.017E+12, 2.018E+12, 1.348E+12, 8.593E+11, 5.236E+11, 3.019E+11, 1.723E+11, 9.614E+10,
 5.272E+10, 2.855E+10, 1.539E+10, 8.098E+09, 4.259E+09, 2.238E+09, 1.187E+09, 6.252E+08,
 3.300E+08, 1.739E+08, 9.118E+07, 4.802E+07, 2.535E+07, 1.338E+07, 7.041E+06, 3.728E+06,
 1.962E+06, 1.039E+06, 5.460E+05, 2.877E+05, 1.518E+05, 7.950E+04, 4.176E+04]


H2 = [1.397E+20, 1.288E+20, 1.186E+20, 1.088E+20, 9.941E+19, 9.234E+19, 8.725E+19, 8.227E+19,
 7.751E+19, 7.286E+19, 6.837E+19, 6.412E+19, 6.000E+19, 5.599E+19, 5.340E+19, 5.089E+19,
 4.852E+19, 4.616E+19, 4.383E+19, 4.164E+19, 3.943E+19, 3.739E+19, 3.529E+19, 3.332E+19,
 3.143E+19, 2.964E+19, 2.763E+19, 2.559E+19, 2.368E+19, 2.187E+19, 2.015E+19, 1.858E+19,
 1.711E+19, 1.578E+19, 1.453E+19, 1.336E+19, 1.223E+19, 1.114E+19, 1.006E+19, 8.977E+18,
 7.929E+18, 6.927E+18, 6.023E+18, 5.230E+18, 3.951E+18, 2.842E+18, 1.988E+18, 1.423E+18,
 1.035E+18, 7.624E+17, 5.396E+17, 3.850E+17, 2.765E+17, 2.000E+17, 1.462E+17, 1.079E+17,
 7.712E+16, 5.547E+16, 4.010E+16, 2.901E+16, 2.103E+16, 1.524E+16, 1.107E+16, 8.030E+15,
 5.835E+15, 4.240E+15, 3.083E+15, 2.242E+15, 1.631E+15, 1.188E+15, 8.647E+14, 6.301E+14,
 4.592E+14, 3.349E+14, 2.445E+14, 1.785E+14, 1.304E+14, 9.551E+13, 6.994E+13, 5.127E+13,
 3.623E+13, 2.567E+13, 1.824E+13, 1.252E+13, 8.329E+12, 5.332E+12, 3.450E+12, 2.243E+12,
 1.468E+12, 9.690E+11, 6.476E+11, 4.320E+11, 2.916E+11, 1.987E+11, 1.373E+11, 9.522E+10,
 6.651E+10, 4.664E+10, 3.277E+10, 2.316E+10, 1.644E+10, 1.170E+10, 8.324E+09, 5.955E+09,
 4.253E+09, 3.051E+09, 2.182E+09, 1.564E+09, 1.123E+09, 8.037E+08, 5.764E+08]



#Juno UVS extended atmospheric profile
data = genfromtxt('../../data/Jupiter_UVS.txt')
H = flip(data[:, 0]) #Km
T0 = flip(data[:, 1]) #K
P0 = flip(data[:, 2]) #barye
Hc = flip(data[:, 3])*1E6 #cc->1/m^3
H2c = flip(data[:, 4])*1E6 #cc->1/m^3
Hec = flip(data[:, 5])*1E6 #cc->1/m^3
CH4c = flip(data[:, 6])*1E6 #cc->1/m^3
C2H2c = flip(data[:, 7])*1E6 #cc->1/m^3
C2H4c = flip(data[:, 8])*1E6 #cc->1/m^3
C2H6c = flip(data[:, 9])*1E6 #cc->1/m^3

OP1 = (a1 + a3)
OP2 = (a2 + a4)
OP5 = a5
OP6 = a6



#Read the wavelength grid from Liu 1995 H2 emission spectrum
UV_dat = genfromtxt('../data/UV_spec_Liu1995.txt')
wav_l = UV_dat[:, 0]*1e-10
RI = UV_dat[:, 1]
RI_max = max(RI)
RI_min = min(RI)
RI_func = interp1d(wav_l, RI, fill_value=(0, 0), bounds_error=False)
RI_norm = RI_func(w_sol)
RI_new = RI_norm/trapz(RI_norm, w_sol/1e-9)

UV_yield = 0.18 #Waite (1983) - Lyman and Werner band partion

O1 = zeros(len(a1))
O2 = zeros(len(a2))
O5 = zeros(len(a5))
O6 = zeros(len(a6))
inx = len(H)
Odepth = zeros([ inx, len(a1)])
UV_emission = zeros([ inx, len(a1)])
inn = where(H <= 150)

fig, ax = plt.subplots(1, 1, figsize=(6, 6))


#CSDA heating rates
#Waite test: Jupiter-Auroral-Energy-Deposition/WaiteHR.txt
m1 = genfromtxt('../../Jupiter-Auroral-Energy-Deposition/Hrates/Benmahi85keV.txt')
H1 = m1[:,0] #Height (km)
N1 = m1[:,1]*1.6e-19*1e6 #Neutral heating rate (eV/cm^3s)-> (W/m^3)
E1 = m1[:,2]*1.6e-19*1e6 #Electron heating rate (eV/cm^3s)-> (W/m^3)
Pgrid = array(flip(Pgrid))
Hgrid = array(flip(Hgrid))
Hfunc = interp1d(Hgrid, Pgrid, fill_value=(0, 0), bounds_error=False)
Nfunc = interp1d(H1, N1, fill_value=(0, 0), bounds_error=False)
#Gerard heating rates
print(H1)

#PLANETOCOSMICS heating rates
m4 = genfromtxt('../../data/PLANETOCOSMICS_noBrhm.txt')
H4 = m4[:, 0]
E4 = m4[:, 1]
g4 = m4[:, 2]
Ev4 = 0.01*E4*1E-3*g4*1e6
P1 = Hfunc(H1)
P4 = Hfunc(H4)
Efunc = interp1d(P4, Ev4, fill_value=(0, 0), bounds_error=False)
Ev5 = Efunc(P1)
N1_new = Nfunc(H)
Tau = zeros(len(w_sol))

for i in range(0, len(H) - 2):
   dz = (H[i] - H[i+1])*1e3
   O1 = (H2c[i]*OP1*dz)
   O2 = (Hec[i]*OP2*dz)
   O5 = (CH4c[i]*OP5*dz)
   O6 = (C2H2c[i]*OP6*dz)
   Tau = Tau + (O1 + O2 + O5 + O6)
   UV_prod1 = RI_new*N1_new[i]*UV_yield*w_sol*1e-4/(0.33*6.626e-34*3e8*1e6)
   UV_prod2 = RI_new*N1_new[i+1]*UV_yield*w_sol*1e-4/(0.33*6.626e-34*3e8*1e6)
   UV_prod = sqrt(UV_prod1*UV_prod2)
   UV_net = UV_prod*exp(-1*Tau)
   UV_net = where(UV_net < 0.0 , 0.0, UV_net)
   UV_emission[i, :] = UV_net
   Odepth[i, :] = Tau

integral_UV = -1*trapz(UV_emission, H*1e3, axis=0)


# Wavelength ranges in nm
range1_min, range1_max = 123, 130  # nm
range2_min, range2_max = 155, 162  # nm

# Convert to meters for comparison
range1_min_m, range1_max_m = range1_min * 1e-9, range1_max * 1e-9
range2_min_m, range2_max_m = range2_min * 1e-9, range2_max * 1e-9

# Find indices for the wavelength ranges
indices_range1 = np.where((w_sol >= range1_min_m) & (w_sol <= range1_max_m))[0]
indices_range2 = np.where((w_sol >= range2_min_m) & (w_sol <= range2_max_m))[0]

# Integrate over the selected wavelength ranges
I_123_130 = np.trapz(integral_UV[indices_range1], w_sol[indices_range1]) / 1e3  # kR
I_155_162 = np.trapz(integral_UV[indices_range2], w_sol[indices_range2]) / 1e3  # kR

# Compute the color ratio
color_ratio = I_155_162 / I_123_130
sum_Intensity = sum(integral_UV)/1e3

print("color ratio: " + str(color_ratio))
print("Total emission (kR) :" + str(trapz(integral_UV, w_sol/1e-9)/1e3 ))

'''
ax.plot(w_sol/1e-9, integral_UV/trapz(integral_UV, w_sol/1e-9), c = 'C1')
ax.plot(w_sol/1e-9, RI_new, c = 'b')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Relative Intensity')
ax.set_yscale('log')
ax.legend(['TOA UV Emission (PJ7)', 'H$_{2}$ emission (Liu et al., 1995)'])
ax.set_xlim([120, 160])
savefig('PJ7_UVEmission.png')
'''






