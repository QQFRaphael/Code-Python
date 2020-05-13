##  This script calculate 2D barotropic Rossby wave
##  ray paths under a specified background flow.
##  Details can be found: Karoly, 1983, Dynamics of Atmospheres and Oceans
##  and Jeff Shaman and Eli Tziperman, 2005, Journal of Climate
## 
##  Written by QQF, @ZJU, 2019-08-06

#############################
##  A) Ks=sqrt(k^2+l^2)
##  Karoly (1983) Eqn 10:
##  B) dk/dt=-k*dUbarM/dx -l*dVbarM/dx +(d2qbar/dxy*k-d2qbar/dxx*l)/Ks^2
##  C) dl/dt=-k*dUbarM/dy -l*dVbarM/dy +(d2qbar/dyy*k-d2qbar/dxy*l)/Ks^2
##  Karoly (1983) Eqn 9:
##  D) dx/dt=ug=UbarM+{(k^2-l^2)*dqbar/dy - 2*k*l*dqbar/dx}/Ks^4
##  E) dy/dt=vg=VbarM+{2*k*l*dqbar/dy - (k^2-l^2)*dqbar/dx}/Ks^4
##  Karoly (1983) Eqn 8:
##  F) omega=UbarM*k+VbarM*l+(dqbar/dx*l-dqbar/dy*k)/Ks^2
##  The above A-F are the key equations to derive wave ray paths
##  I use Runge-Kutta method to perform integration...

from namelist import *
from util_func import *
from scipy import interpolate


datapath = "./data/input.nc"
q, um, vm, umx, vmx, umy, vmy, qx, qy, qxx, qxy, qyy, lat, lon, xm, ym, xm360 = read(datapath)

nlat = np.shape(lat)[0]
nlon = np.shape(lon)[0]

#print(q, um, vm, umx, vmx, umy, vmy, qx, qy, qxx, qxy, qyy, lat, lon, xm, ym, xm360)

# calculate some derivatives
um_int_func = interpolate.interp2d(xm, ym, um, kind='linear')
vm_int_func = interpolate.interp2d(xm, ym, vm, kind='linear')
umx_int_func = interpolate.interp2d(xm, ym, umx, kind='linear')
umy_int_func = interpolate.interp2d(xm, ym, umy, kind='linear')
vmx_int_func = interpolate.interp2d(xm, ym, vmx, kind='linear')
vmy_int_func = interpolate.interp2d(xm, ym, vmy, kind='linear')
qx_int_func = interpolate.interp2d(xm, ym, qx, kind='linear')
qy_int_func = interpolate.interp2d(xm, ym, qy, kind='linear')
qxx_int_func = interpolate.interp2d(xm, ym, qxx, kind='linear')
qxy_int_func = interpolate.interp2d(xm, ym, qxy, kind='linear')
qyy_int_func = interpolate.interp2d(xm, ym, qyy, kind='linear')



spotk = k_wavenumbers / r / np.cos(frcy)

##  Calculate the initial l wave number from the initial omega
##  and k by solving the polynomial equation based on the
##  dispersion relation (Eqn 8 in Karoly 1983):
##  hange the following to have a non zero frequency:
coeff = np.zeros(4)
coeff[0] = vm_int_func(frcxm,frcym)[0]
coeff[1] = um_int_func(frcxm,frcym)[0]*spotk-frequency
coeff[2] = vm_int_func(frcxm,frcym)[0]*spotk*spotk+qx_int_func(frcxm,frcym)[0]
coeff[3] = um_int_func(frcxm,frcym)[0]*spotk*spotk*spotk-qy_int_func(frcxm,frcym)[0]*spotk-frequency*spotk*spotk

lroot = np.roots(coeff)

print("initial k = ", spotk)
print("initial l = ", lroot)

# above is initialization, now we start ray tracing

# loop over all l
for idxl in range(0,3):

	spotl = lroot[idxl]
	print("Root # ", idxl, "  spotl = ", spotl)

	if np.not_equal(np.imag(spotl),0): continue

	## Starting the loop with the above initial k,l, and Ks
	lonn = np.empty(Nsteps+1)
	latn = np.empty(Nsteps+1)
	xn = np.empty(Nsteps+1)
	yn = np.empty(Nsteps+1)
	kn = np.empty(Nsteps+1)
	ln = np.empty(Nsteps+1)

	lonn[:] = np.nan
	latn[:] = np.nan
	xn[:] = np.nan
	yn[:] = np.nan
	kn[:] = np.nan
	ln[:] = np.nan

	# start integration
	for t in range(0,Nsteps):
		if np.equal(np.remainder(t,10),0):
			print("step = # ", t)

		if t == 0 :
			x0 = xn[0] = frcxm
			y0 = yn[0] = frcym
			k0 = kn[0] = spotk
			l0 = ln[0] = np.real(spotl)
			lonn[0] = x2lon(frcxm)
			latn[0] = y2lat(frcym)
		else:
			x0, y0, k0, l0 = xn[t], yn[t], kn[t], ln[t]

		dx, dy, dk, dl = RK(x0,y0,k0,l0, \
			um_int_func,vm_int_func, \
			umx_int_func,umy_int_func, \
			vmx_int_func,vmy_int_func, \
			qx_int_func,qy_int_func, \
			qxy_int_func,qxx_int_func,qyy_int_func)

		tn = t+1
		xn[tn] = x0 + dx
		if xn[tn] >= xm360: xn[tn] = xn[tn] - xm360
		yn[tn] = y0 + dy
		kn[tn] = k0 + dk
		ln[tn] = l0 + dl

		# convert Mercator location to lat-lon location
		lonn[tn] = x2lon(xn[tn])
		latn[tn] = y2lat(yn[tn])

	# write out all the results
	if frequency == 0:
		outfile = "./output/raypath_lat_{:.2f}_lon_{:.2f}_period_{}_k_{:d}_root_{:d}.txt".format(frcy/np.pi*180.0,frcx/np.pi*180.0,'inf',k_wavenumbers,idxl)
	else:
		outfile = "./output/raypath_lat_{:.2f}_lon_{:.2f}_period_{:.2f}_k_{:d}_root_{:d}.txt".format(frcy/np.pi*180.0,frcx/np.pi*180.0,2.0*np.pi/(frequency*day),k_wavenumbers,idxl)

	fout = open(outfile,"w")
	#fout.write("latitude        longitude \n")
	fmt = "{:.4f}        {:.4f} \n"
	for idx in range(0,Nsteps):
		fout.write(fmt.format(latn[idx], lonn[idx]))
	fout.close()
