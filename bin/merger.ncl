begin
;--------------------------------------
; READ INPUT
; -------------------------------------
;system("sort -n -k2 output/raincell/irt_tracks_output_pure.txt > output/raincell/irt_tracks_output_pure_sort.txt")
obj = readAsciiTable("output/raincell/irt_objects_output.txt",23 ,"integer",0)
track = readAsciiTable("output/raincell/irt_tracks_output_pure_sort.txt",12 , "integer",0)
track2 = readAsciiTable("input/cp/irt_tracks_header_sorted.txt",5,"integer",0)

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
outvar = (/ispan(1,na,1),ispan(1,na,1),ispan(1,na,1)/)
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

; max size
do i = 0, na-1
  maxarea = max(track(ind(track(:,0) .eq. ta(i)),4))
if (i .gt. 0)
;print(dimsizes(tmaxarea))
;print(dimsizes(track(ind(track(:,0) .eq. ta(i) .and. track(:,4) .eq. maxarea  ),1)))
end if
  tmaxaread = track(ind(track(:,0) .eq. ta(i) .and. track(:,4) .eq. maxarea  ),1)
  if (i .gt. 1)
    tmaxarea = min(track(ind(track(:,0) .eq. ta(i) ),1));just first ;tmaxaread(0)
;print("greater : "+ dimsizes(tmaxarea))
  else 
   tmaxarea = track(ind(track(:,0) .eq. ta(i)),1)
  end if
  delete(tmaxaread)
  
  outvar(2,ta(i)-1) = tmaxarea(0) 

end do

 matrix =    str_insert( outvar(0,:)," ",7) + \
            str_insert(outvar(1,:)," ",7) ;+ \
       ;     str_insert(outvar(2,:)," ",7) 
        n = dimsizes(matrix)
        n(0) = n(0)+1
        print(n)
        outputtxt = new(n,"string")
        outputtxt(0)    = " "
        outputtxt(1::)  = matrix

        asciiwrite("mergingCPs.txt",matrix);zu pdf zusammenfügen


end
