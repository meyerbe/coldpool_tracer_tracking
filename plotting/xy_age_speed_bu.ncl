begin
;--------------------------------
;--- GET DATA
;------------------------------
data = readAsciiTable("../output/cp/coldpool_tracer_out.txt",12,"float",0)

starttime = data(0,0)
ntime = count_unique_values(data(:,0))   ; number of timesteps
n_CPs = count_unique_values(data(:,3))   ; number of CPs for outputdimension
CPs = get_unique_values(data(:,3))       ; CP IDs
maxage = tointeger(max(data(:,1)))       ; oldest CP for output dimension
;--------------------------------
;--- ALLOCATION and INITIALIZATION
;------------------------------

outmean = new((/n_CPs,maxage/),float)
outstd = new((/n_CPs,maxage/),float)

do i = 0,20; ntime -1
  indextime := ind(data(:,0) .eq. starttime+i)
  current_data := data(indextime,:)
  current_CPs := get_unique_values(current_data(:,3))
  c_nCPs = count_unique_values(current_data(:,3))
  do j = 0,c_nCP-1
    cpindex = ind(CPs .eq.current_CPs(j))
    age := current_data(min(cpindex),1) 
    FF := sqrt(current_data(cpindex,10)^2 + current_data(cpindex,11)^2) 
    outmean(cpindex,age-1) = avg(FF)
    outstd(cpindex,age-1) = stddev(FF)
  end do
end do
  
  
end


