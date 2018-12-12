begin
;--------------------------------------
; READ INPUT
; -------------------------------------
;system("sort -n -k2 output/raincell/irt_tracks_output_pure.txt > output/raincell/irt_tracks_output_pure_sort.txt")
path=path
tr=tr ; trashold of precipitation size to start tracer
dx=dx ; domainsize
dy=dy
print(tr)
print(dx)
obj = readAsciiTable(path+"output/raincell/irt_objects_output.txt",23 ,"integer",0)
track = readAsciiTable(path+"output/raincell/irt_tracks_output_pure_sort.txt",12 , "integer",0)
track2 = readAsciiTable(path+"input/cp/irt_tracks_header_sorted.txt",5,"integer",0)

; -----------------------------
; sort data 
; -----------------------------
;tracks terminate as merger 
tIDmall = track2(ind((track2(:,4) .eq. 1)),0)
tIDm = get_unique_values(tIDmall)
ni = count_unique_values(tIDm)
tIDfall = track2(ind((track2(:,3) .eq. 1)),0)
tIDf = get_unique_values(tIDfall)
nj = count_unique_values(tIDf)   
na = count_unique_values(track2(:,0)) ; all events 
ta = get_unique_values(track2(:,0)) 
;-------------------------
; find corresponding IDs
;-------------------------
outvar = (/ispan(1,na,1),ispan(1,na,1),ispan(1,na,1),ispan(1,na,1)/)
; falsch outvar = conform_dims((/na,4/),ispan(1,na,1),0) 
do i=0,ni-1
  ;get the last object of this track, which is the one which terminates by contributing to a larger one
  oID2m=track(max(ind((track(:,0) .eq. tIDm(i)))),3) ;get the last object of this track, which is the one which terminates by contributing to a larger one
  ; get into which object this object merges at the next timestep
  oID2i=obj(ind(obj(:,2) .eq. oID2m),15)
  ; get the corresponding track to this ID, track ID in

  tIDin = track(ind(track(:,3) .eq. oID2i),0)
  outvar(0:1,tIDm(i)-1) = (/tIDm(i),tIDin/)
end do
do i=0,nj-1
  outvar(0:1,tIDf(i)-1) = (/tIDf(i),0/)
print("fragm " + tIDf(i))
end do

;----------------------------------------------------
; identify timestep when tracer should start
; either when gpsize> 50 or at time of maximum size 
;----------------------------------------------------

do i = 0, na-1  ; na is number of all tracks
  ; get the size at every timestep for the current ID i 
  currentareas := track(ind(track(:,0) .eq. ta(i)),4)
  ; first check if precip event ever reaches size above given trashold ... 
  if any(currentareas .gt. tr)  ; and find first timestep when area is above tr
    start = min(track(ind(track(:,0) .eq. ta(i) .and. track(:,4) .gt. tr),1))
  else ; or find timestep when area reaches its maximum
    startarea := max(track(ind(track(:,0) .eq. ta(i)),4))
    start = min(track(ind(track(:,0) .eq. ta(i) .and. track(:,4) .eq. startarea),1))
  end if
  outvar(2,ta(i)-1) = start 
  outvar(3,ta(i)-1) = min(track(ind(track(:,0) .eq. ta(i)),1))
end do

; -------------------------------------------------
; DATA OUTPUT 
; -----------------------------------------------
 matrix =    str_insert( outvar(0,:)," ",7) + \
            str_insert(outvar(1,:)," ",7) + \
            str_insert(outvar(2,:)," ",7)+ \
            str_insert(outvar(3,:)," ",7)

        n = dimsizes(matrix)
        n(0) = n(0)+1
        print(n)
        outputtxt = new(n,"string")
        outputtxt(0)    = " "
        outputtxt(1::)  = matrix

        asciiwrite(path+"input/cp/mergingCPs.txt",matrix);zu pdf zusammenfügen

system("echo "+max(ta) +" > na.txt")

end
