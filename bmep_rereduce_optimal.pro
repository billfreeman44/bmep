;program to rereduce LRIS data if the data was saved in a SAVE FILE!!
;this method is largely obsolete.
;working
pro bmep_rereduce_optimal,savefile
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
    
  ;restore the file
  if file_test(savefile) then restore,savefile else message,'no save file'
  
  ;if it has the state structure then rereduce.
  if n_elements(state) ne 0 then begin
    ;extraction!!
    bmep_extraction,j,centerarr,state,order,widtharr,bkgndl,bkgndr,$
      printp,fitgaussp,plotp,slidep,pwindowsize,singlep,$
      cosmic_sigma,n_iterate_cosmic,bkgnd_naverage,max_rays_per_col,$
      $;OUTPUTS
      F,ferr,Fopt,fopterr,img2d_nobkgnd,parr,sky_reslut_arr,sky_residuals,revisevar=state.revisevar
      
    ;do not plot bad pixels
    bmep_calc_skyline_mask_v2,ferr,skymaskpercent,skymask
    
    ;individual post reduction processing!!
    if state.raw_slitname eq 'mips313' then begin
      print,'editing mips313'
      index=where(wavel eq 1736.16,ct)
      if ct ne 0 then begin
        print,'edit successful'
        F[index]=0.0
        Fopt[index]=0.0
        ferr[index]=0
        fopterr[index]=0
      endif
    endif
    
    if state.raw_slitname eq 'mips358' then begin
      print,'editing mips358'
      index=where(wavel eq 2639.50,ct)
      if ct ne 0 then  begin
        print,'edit successful'
        F[index]=0.0
        Fopt[index]=0.0
        ferr[index]=0
        fopterr[index]=0
      endif
    endif
    
    if state.raw_slitname eq 'mips374' then begin
      print,'editing mips374'
      index=where(wavel eq 2768.55,ct)
      if ct ne 0 then begin
        print,'edit successful'
        F[index]=0.0
        Fopt[index]=0.0
        ferr[index]=0
        fopterr[index]=0
      endif
    endif
    
    
    
    flux=f
    err=ferr
    flux_opt=fopt
    erropt=fopterr
    wavel=state.wavel
    slitname=state.slitname
    img2d=state.data
    flux_opt_masked=fopt*skymask
    flux_masked=f*skymask
    
    ;save the newly extracted spectra
    if n_elements(centerarr) eq 1 then suffix='' else suffix='---'+ssi(j+1)
    save,/variables,filename=savefile
    print,'saved ',slitname+'.sav'
    
    ;save the newly extracted spectra in text format
    forprint,wavel,f,ferr,textout=strmid(savefile,0,strlen(savefile)-4)+'.txt',/nocomment
    forprint,wavel,Fopt,Fopterr,textout=strmid(savefile,0,strlen(savefile)-4)+'_optimal.txt',/nocomment
    
  endif else begin
    print,'COULDNT rereduc ',savefile
  endelse
  print,'done rereducing ',savefile
end

