begin
begTime = get_cpu_time()
; -------------- SETTINGS -----------
back = 10 ; search bakwards for max updraft by "back" timesteps

; -------------- READ VARIABLES ---------------------

mask_file =  addfile("../output/raincell/irt_tracks_mask.nc", "r")
w_file =  addfile("/nbi/ac/conv1/henneb/Moseley2018/data/lind_p2K/lind_p2K.out.vol.w.nc","r")

rmask = mask_file->var1(:,0,:,:)
wmask = w_file->w({rmask&time(tstartind)-back*dt+1:rmask&time(tstartind) },:,:,:)

do t = 0,5
 
  rmask = mask_file->var1(t,0,:,:)
  rtimestep = rmask&time*24*60/5. + 1
  wtemp = w_file->w(rtimestep-back:rtimestep,:,:,{2000-3000})
  





end do





