import numpy as np
from namelist import *
from netCDF4 import Dataset


def read(filepath):
	nc = Dataset(filepath)
	q = np.array(nc.variables['q'][:])
	um = np.array(nc.variables['um'][:])
	vm = np.array(nc.variables['vm'][:])
	umx = np.array(nc.variables['umx'][:])
	vmx = np.array(nc.variables['vmx'][:])
	umy = np.array(nc.variables['umy'][:])
	vmy = np.array(nc.variables['vmy'][:])
	qx = np.array(nc.variables['qx'][:])
	qy = np.array(nc.variables['qy'][:])
	qxx = np.array(nc.variables['qxx'][:])
	qxy = np.array(nc.variables['qxy'][:])
	qyy = np.array(nc.variables['qyy'][:])
	lat = np.array(nc.variables['lat'][:])
	lon = np.array(nc.variables['lon'][:])
	xm = np.array(nc.variables['xm'][:])
	ym = np.array(nc.variables['ym'][:])
	xm360 = np.array(nc.variables['xm360'][:])

	return q, um, vm, umx, vmx, umy, vmy, qx, qy, qxx, qxy, qyy, lat, lon, xm, ym, xm360


def ug(k, l, um, qx, qy):
	"""
	This is Eqn.9, ug in Karoly, 1983
	"""
	K2 = k*k + l*l
	K4 = K2 * K2
	return um + ((k*k-l*l)*qy - 2.0*k*l*qx)/K4


def vg(k, l, vm, qx, qy):
	"""
	This is Eqn.9, vg in Karoly, 1983
	"""
	K2 = k*k + l*l
	K4 = K2 * K2
	return vm + (2.0*k*l*qy + (k*k-l*l)*qx)/K4


def dk_dt(k, l, umx, vmx, qxy, qxx):
	"""
	This is Eqn.10, dk/dt in Karoly, 1983
	"""
	K2 = k*k + l*l
	K4 = K2 * K2
	return -k*umx - l*vmx + (qxy*k-qxx*l)/K4


def dl_dt(k, l, umy, vmy, qxy, qyy):
	"""
	This is Eqn.10, dl/dt in Karoly, 1983
	"""
	K2 = k*k + l*l
	K4 = K2 * K2
	return -k*umy - l*vmy + (qyy*k-qxy*l)/K4


def y2lat(yy):
	"""
	Mercator y to lat
	"""
	return 180.0/np.pi*(2.0*np.arctan(np.exp(yy/r))-np.pi/2.0)


def x2lon(xx):
	"""
	Mercator x to lon
	"""
	return xx*180.0/np.pi/r


def RK(x0,y0,k0,l0,um_int_func,vm_int_func,umx_int_func,umy_int_func,vmx_int_func,vmy_int_func,qx_int_func,qy_int_func,qxy_int_func,qxx_int_func,qyy_int_func):
	"""
	This function perform Runge-Kutta integration
	"""
	# RK step1
	kx0 = ug(k0,l0,um_int_func(x0,y0)[0],qx_int_func(x0,y0)[0],qy_int_func(x0,y0)[0])
	ky0 = vg(k0,l0,vm_int_func(x0,y0)[0],qx_int_func(x0,y0)[0],qy_int_func(x0,y0)[0])
	kk0 = dk_dt(k0,l0,umx_int_func(x0,y0)[0],vmx_int_func(x0,y0)[0],qxy_int_func(x0,y0)[0],qxx_int_func(x0,y0)[0])
	kl0 = dl_dt(k0,l0,umy_int_func(x0,y0)[0],vmy_int_func(x0,y0)[0],qxy_int_func(x0,y0)[0],qyy_int_func(x0,y0)[0])

	# RK step2
	x1 = x0 + 0.5*kx0*dt
	y1 = y0 + 0.5*ky0*dt
	k1 = k0 + 0.5*kk0*dt
	l1 = l0 + 0.5*kl0*dt

	kx1 = ug(k1,l1,um_int_func(x1,y1)[0],qx_int_func(x1,y1)[0],qy_int_func(x1,y1)[0])
	ky1 = vg(k1,l1,vm_int_func(x1,y1)[0],qx_int_func(x1,y1)[0],qy_int_func(x1,y1)[0])
	kk1 = dk_dt(k1,l1,umx_int_func(x1,y1)[0],vmx_int_func(x1,y1)[0],qxy_int_func(x1,y1)[0],qxx_int_func(x1,y1)[0])
	kl1 = dl_dt(k1,l1,umy_int_func(x1,y1)[0],vmy_int_func(x1,y1)[0],qxy_int_func(x1,y1)[0],qyy_int_func(x1,y1)[0])

	# RK step3
	x2 = x0 + 0.5*kx1*dt
	y2 = y0 + 0.5*ky1*dt
	k2 = k0 + 0.5*kk1*dt
	l2 = l0 + 0.5*kl1*dt

	kx2 = ug(k2,l2,um_int_func(x2,y2)[0],qx_int_func(x2,y2)[0],qy_int_func(x2,y2)[0])
	ky2 = vg(k2,l2,vm_int_func(x2,y2)[0],qx_int_func(x2,y2)[0],qy_int_func(x2,y2)[0])
	kk2 = dk_dt(k2,l2,umx_int_func(x2,y2)[0],vmx_int_func(x2,y2)[0],qxy_int_func(x2,y2)[0],qxx_int_func(x2,y2)[0])
	kl2 = dl_dt(k2,l2,umy_int_func(x2,y2)[0],vmy_int_func(x2,y2)[0],qxy_int_func(x2,y2)[0],qyy_int_func(x2,y2)[0])

	# RK step4
	x3 = x0 + kx2*dt
	y3 = y0 + ky2*dt
	k3 = k0 + kk2*dt
	l3 = l0 + kl2*dt

	kx3 = ug(k3,l3,um_int_func(x3,y3)[0],qx_int_func(x3,y3)[0],qy_int_func(x3,y3)[0])
	ky3 = vg(k3,l3,vm_int_func(x3,y3)[0],qx_int_func(x3,y3)[0],qy_int_func(x3,y3)[0])
	kk3 = dk_dt(k3,l3,umx_int_func(x3,y3)[0],vmx_int_func(x3,y3)[0],qxy_int_func(x3,y3)[0],qxx_int_func(x3,y3)[0])
	kl3 = dl_dt(k3,l3,umy_int_func(x3,y3)[0],vmy_int_func(x3,y3)[0],qxy_int_func(x3,y3)[0],qyy_int_func(x3,y3)[0])

	#RK results
	dx = dt * (kx0+2*kx1+2*kx2+kx3)/6
	dy = dt * (ky0+2*ky1+2*ky2+ky3)/6
	dk = dt * (kk0+2*kk1+2*kk2+kk3)/6
	dl = dt * (kl0+2*kl1+2*kl2+kl3)/6

	return dx, dy, dk, dl