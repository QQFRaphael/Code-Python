load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

; NCEP DOE Reanalysis II data
; uwnd, vwnd are input, and the output is stream function
; Note 90N and 90S should be eliminated for the convenience of following calculation

system("rm -rf input.nc")

; JFM
idxA = 0
idxB = 1
idxC = 2

mylev = 200
yearstart = 1979
yearend = 2018 
nyear = yearend - yearstart + 1

ufile = addfile("uwnd.mon.mean.nc", "r")
vfile = addfile("vwnd.mon.mean.nc", "r")

lev = ufile->level
lev_idx = ind(lev.eq.mylev)

lat = ufile->lat
lon = ufile->lon

nlat = dimsizes(lat)
nlon = dimsizes(lon)

year=calendar_decode2(ufile->time,0)
year_idx=ind(year(:,0).ge.(yearstart).and.year(:,0).le.(yearend))  

uwnds = short2flt(ufile->uwnd(year_idx,lev_idx,:,:))
vwnds = short2flt(vfile->vwnd(year_idx,lev_idx,:,:))

uA = uwnds(idxA::12,:,:)
uB = uwnds(idxB::12,:,:)
uC = uwnds(idxC::12,:,:)
u = (uA + uB + uC) / 3.0
copy_VarMeta(uB, u)

vA = vwnds(idxA::12,:,:)
vB = vwnds(idxB::12,:,:)
vC = vwnds(idxC::12,:,:)
v = (vA + vB + vC) / 3.0
copy_VarMeta(vB, v)

uwnd = dim_avg_n_Wrap(u,0)
vwnd = dim_avg_n_Wrap(v,0) 

pi = atan(1.0)*4. 
re=6378388 

vort = uv2vr_cfd(uwnd,vwnd,lat,lon,3)
q = vort + conform(vort,2.0*7.2925e-5*sin(lat*pi/180.0),0)
copy_VarMeta(uwnd, q)

lat_r = lat * pi / 180.0
lon_r = lon * pi / 180.0

xm = re * lon_r 
ym = re * log( (1+sin(lat_r)) / cos(lat_r))

xm360 = 360.0 * re * pi / 180.0

coslat = cos(lat_r)
cos2 = coslat * coslat

u_xy = grad_latlon_cfd(uwnd, lat, lon, True, False)
uy = u_xy[0]
ux = u_xy[1]

v_xy = grad_latlon_cfd(vwnd, lat, lon, True, False)
vy = v_xy[0]
vx = v_xy[1]

q_xy = grad_latlon_cfd(q, lat, lon, True, False)
qy = q_xy[0]
qx = q_xy[1]

qq_xy = grad_latlon_cfd(qx, lat, lon, True, False)
qxy = q_xy[0]
qxx = q_xy[1]

qqq_xy = grad_latlon_cfd(qy, lat, lon, True, False)
qyy = q_xy[0]
qyx = q_xy[1]

if(any(ismissing(uy))) then
    print("uy contains some missing values. Beware.")
end if

if(any(ismissing(ux))) then
    print("ux contains some missing values. Beware.")
end if

if(any(ismissing(vy))) then
    print("vy contains some missing values. Beware.")
end if

if(any(ismissing(vx))) then
    print("vx contains some missing values. Beware.")
end if

if(any(ismissing(qy))) then
    print("qy contains some missing values. Beware.")
end if

if(any(ismissing(qx))) then
    print("qx contains some missing values. Beware.")
end if

if(any(ismissing(qxy))) then
    print("qxy contains some missing values. Beware.")
end if

if(any(ismissing(qxx))) then
    print("qxx contains some missing values. Beware.")
end if

if(any(ismissing(qyy))) then
    print("qyy contains some missing values. Beware.")
end if

um = uwnd / conform(uwnd, cos(lat_r), 0)
vm = vwnd / conform(uwnd, cos(lat_r), 0)
q = q / conform(uwnd, cos(lat_r), 0)

umx = ux / conform(uwnd, cos(lat_r), 0)
vmx = vx / conform(uwnd, cos(lat_r), 0)

umy = uy / conform(uwnd, cos(lat_r), 0)
vmy = vy / conform(uwnd, cos(lat_r), 0)

qx = qx / conform(uwnd, cos(lat_r), 0)
qy = qy / conform(uwnd, cos(lat_r), 0)

qxx = qxx / conform(uwnd, cos(lat_r), 0)
qxy = qxy / conform(uwnd, cos(lat_r), 0)
qyy = qyy / conform(uwnd, cos(lat_r), 0)

copy_VarMeta(uwnd, q)
copy_VarMeta(uwnd, um)
copy_VarMeta(uwnd, vm)
copy_VarMeta(uwnd, umx)
copy_VarMeta(uwnd, vmx)
copy_VarMeta(uwnd, umy)
copy_VarMeta(uwnd, vmy)
copy_VarMeta(uwnd, qx)
copy_VarMeta(uwnd, qy)
copy_VarMeta(uwnd, qxy)
copy_VarMeta(uwnd, qxx)
copy_VarMeta(uwnd, qyy)
copy_VarMeta(lat, ym)
copy_VarMeta(lon, xm)

; do not include 90S and 90N for the convenient of computation
out = addfile("input.nc","c")
out->lat = lat(1:nlat-2) 
out->lon = lon
out->q = q(1:nlat-2,:)
out->um = um(1:nlat-2,:)
out->vm = vm(1:nlat-2,:)
out->umx = umx(1:nlat-2,:)
out->umy = umy(1:nlat-2,:)
out->vmx = vmx(1:nlat-2,:)
out->vmy = vmy(1:nlat-2,:)
out->qx = qx(1:nlat-2,:)
out->qy = qy(1:nlat-2,:)
out->qxx = qxx(1:nlat-2,:)
out->qxy = qxy(1:nlat-2,:)
out->qyy = qyy(1:nlat-2,:)
out->xm = xm
out->ym = ym(1:nlat-2)
out->xm360 = xm360
