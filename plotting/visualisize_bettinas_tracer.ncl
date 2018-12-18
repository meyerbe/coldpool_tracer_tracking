begin
; SETTINGS 
startzeit =1 ; time when first CP is tracked 

;--------------------------------
;--- GET DATA
;------------------------------
; get tracer heade to identify self imerged events
; expid ="single_3D_noise/dTh3_z1000_r1000"
 expid =expid ;"3D_sfc_fluxes_off/triple_3D_noise/Out_CPDry_triple_dTh3K"
 lvl=lvl
;--- tracer data
 data_all   = readAsciiTable(expid+"/tracer_k"+lvl+"/output/coldpool_tracer_out.txt",18,"float",0) 
 print(expid+"/tracer_k"+lvl+"/output/coldpool_tracer_out.txt")

 timei=0
 agei=1
 traceIDi= 2
 cpIDi = 3
 Xpos = 4
 Ypos = 5
 Xgp = 6
 Ygp = 7
 data = data_all;           (ind(data_all(:,cpIDi) .eq. any(track0)),:)
 print("number of tracers: " +max(data(:,traceIDi)))
 print("number of cps: "+count_unique_values(data_all(:,cpIDi)))
 print("number of timesteps: " +max(data(:,timei)))
 print("max CP age: " +max(data(:,agei)))

 print("dims")
;-------------------
;---modelo output
 filein = addfile(expid+"/tracer_k"+lvl+"/input/uv_alltimes.nc","r")
 dna = getfilevardims(filein,"u")
 ;u = filein->u($dna(0)$ | :, $dna(2)$  | :, $dna(1)$ | :) 
 ;v = filein->v($dna(0)$ | :, $dna(2)$  | :, $dna(1)$ | :)
 u = filein->u( :, :, :) 
 v = filein->v( :, :, :)
 print(u!0)
 print(u!1)
 print(u!2)
 dims = dimsizes(u)
 resv =True
 resv@gsnFrame     = False
 resv@gsnDraw      = False
;  resv@gsnMaximize  = True     ; Maximize plot in frame.
;  resv@trGridType = "TriangularMesh"  ;makes it faster
;  resv@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last
;  resv@vpClipOn     = False    ; vector over plot border
  resv@vcMinDistanceF             = 0.01
  resv@vcRefLengthF               = 0.05
  resv@vcRefMagnitudeF         = 5.0              ; make vectors larger
  resv@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors;   
  resv@vcPositionMode ="ArrowHead" ;"ArrowCenter"; "ArrowHead";Default: ArrowCenter (/"ArrowCenter","ArrowTail","ArrowHead"/)
  resv@tiMainOn = False
  ;resv@tiMainString = "Vector plot with default clipping"
  resv@vcRefAnnoOn = False
  resv@vcRefAnnoJust = 0
  resv@vcRefAnnoSide ="Right"
  resv@vcRefAnnoParallelPosF = 1-0.6
  resv@vcRefAnnoOrthogonalPosF = -0.3
 ; resv@trXMaxF =350 ;250
 ; resv@trYMaxF =350; ;300;
 ; resv@trXMinF =50;220
 ; resv@trYMinF =50;270;


;--- select CP Tracer
;find first CPid:
 cpIDs = get_unique_values(data(:,cpIDi))      ; individual CP numbers
print("printing cpIDs: "+cpIDs) 
noCPs = count_unique_values(data(:,cpIDi))    ; number od CPs 
 notracer = count_unique_values(data(:,traceIDi))

 colors= (/"orange","darkgreen","magenta","blue","red","seagreen"/)
 colortables=(/"MPL_Blues","MPL_Greens","MPL_Purples","MPL_Oranges","MPL_Reds","MPL_Greys","MPL_PuRd"/)
 marker=new((/noCPs/),graphic)
 markers=new((/noCPs/),graphic)
 mres = True
 mres@gsEdgeColor  = "Black"
 cmap =read_colormap_file("default");uniform");BlAqGrYeOrReVi200")

 print("start main loop")
 do t = 0,dims(0)-1,1 
   wks = gsn_open_wks("png",expid+"/tracer_k"+lvl+"/plot_"+sprinti("%0.4i",t) )
   ;pl = gsn_csm_vector(wks, u(t,  :,  :),v(t, :,:), resv)
   pl = gsn_csm_vector(wks, u($dna(0)$|t, $dna(2)$|:, $dna(1)$|:),v($dna(0)$|t, $dna(2)$|:, $dna(1)$|:), resv)

   if (.not. all(ismissing(ind(data(:,0) .eq. t+startzeit))))
     datatime = data(ind(data(:,0) .eq. t+startzeit),:) 
     c = 0
     do i = 0,noCPs-1 ; loop troughy cps 
      indcpID = ind(datatime(:,cpIDi) .eq. cpIDs(i))
      ; print(indcpID)
      mres@gsMarkerColor =  cmap(mod(i,30),:)
      xpt =  datatime(indcpID,Xpos)
      ypt=  datatime(indcpID,Ypos)
      mres@gsMarkerIndex = 16;j
      mres@gsMarkerSizeF = 4
      marker(i) = gsn_add_polymarker(wks,pl,xpt,ypt,mres) 
      delete([/xpt,ypt/])

      mkres               = True
      mkres@gsMarkerIndex = 16
      mres@gsMarkerSizeF = 10 
      mkres@gsMarkerColor = cmap(mod(i,30),:)
      ;markers(i) = gsn_add_polymarker(wks,pl,cogx,cogy,mkres)
      mkres@txFontColor = cmap(mod(i,30),:)
      mkres@txFontHeightF =0.02 
      ;mkres@txFontThicknessF = 0.5
      delete([/indcpID/])
     end do
     delete(datatime)
     draw(pl)
     frame(wks)
   end if
 end do
end
