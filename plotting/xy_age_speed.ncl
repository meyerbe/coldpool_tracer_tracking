begin
;--------------------------------
;--- GET DATA
;------------------------------
data = readAsciiTable("../output/cp/coldpool_tracer_out.txt",12,"float",(/0,5000/))

starttime = data(0,0)
ntime =112; count_unique_values(data(:,0))   ; number of timesteps
n_CPs = 200;count_unique_values(data(:,3))   ; number of CPs for outputdimension
CPs = get_unique_values(data(:,3))       ; CP IDs
maxage = tointeger(max(data(:,1)))       ; oldest CP for output dimension
;--------------------------------
;--- ALLOCATION and INITIALIZATION
;------------------------------
outmean = new((/n_CPs,maxage/),float)
outstd = new((/n_CPs,maxage/),float)

do i = 0,100; ntime -1
  indextime := ind(data(:,0) .eq. starttime+i)
  current_data := data(indextime,:)
  current_CPs := get_unique_values(current_data(:,3))
  c_nCPs = count_unique_values(current_data(:,3))
  do j = 0,min((/c_nCPs-1,6/))
    cpcount := ind(CPs .eq.current_CPs(j))
    cpindex := ind(current_data(:,3) .eq. current_CPs(j))
    age := tointeger(current_data(min(cpindex),1) )
    FF := sqrt(current_data(cpindex,10)^2 + current_data(cpindex,11)^2) 
    outmean(cpcount,age-1) = avg(FF)
    outstd(cpcount,age-1) = stddev(FF)
  end do
end do
  
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("png","gsn_xy")                       ; send graphics to PNG file   
  cmap       = read_colormap_file("srip_reanalysis")
  res                        = True                        ; plot mods desired
  res@vpHeightF     = 0.4
  res@vpWidthF     = 0.6
  res@xyLineColors  = cmap(::2,:);(/"blue","sienna","forestgreen","hotpink","orange",/) 
  res@xyLineThicknesses = conform_dims(10,4,0)
  res@xyDashPatterns   = conform_dims(10,0,0)
  do i =0,10
    plot=gsn_y(wks,outmean(10*i:10*i+9,0:15),res)
;    plot=gsn_y(wks,outstd(10*i:10*i+9,0:15),res)

  end do
  
end


