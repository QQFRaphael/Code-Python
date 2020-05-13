load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

npoints = 360

filename = "raypath_lat_45.00_lon_210.00_period_inf_k_4_root_0"

locs = asciiread("./output/"+filename+".txt", (/npoints,2/), "float")
lat = locs(:,0)
lon = locs(:,1)

res = True
res@gsnDraw=False
res@gsnFrame=False
res@gsnMaximize=True


;res@mpLimitMode = "LatLon"
;res@mpMinLatF = 10
;res@mpMaxLatF = 80
;res@mpMinLonF = 180
;res@mpMaxLonF = 360


;res@mpCenterLonF=180

res@mpFillOn=False

res@lbOrientation="Vertical"

res@gsnLeftString=""
res@gsnRightString=""

wks = gsn_open_wks("png", filename)

map = gsn_csm_map_ce(wks, res)

resc=True
resc@gsLineColor      = "red"                     ; color of lines
resc@gsLineThicknessF = 6.0                       ; thickness of lines

dumm = new(npoints,graphic)
do i = 0 , npoints-2
    dumm(i)=gsn_add_polyline(wks,map,lon(i:i+1),lat(i:i+1),resc)
end do

draw(map)
frame(wks)
