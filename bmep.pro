;version 1.00
;;git commit -a -m ''
;extract data multiple times using different widths to find the best SNR
pro bmep_auto_width_calculator,j,centerarr,state,order,bkgndl,bkgndr,$
    printp,fitgaussp,plotp,slidep,pwindowsize,singlep,cosmic_sigma,n_iterate_cosmic,bkgnd_naverage,max_rays_per_col,$
    $;OUTPUTS
    F,ferr,Fopt,fopterr,img2d_nobkgnd,parr,sky_reslut_arr,sky_residuals,$
    widtharr,leftstatspos,rightstatspos
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  print,'\begin{tabular}{| l | l | l | l | l | l |}'
  print,'  \hline'
  print,'Width & SNR boxcar & SNR optimal   & mean box & mean opt & total flux \\ \hline'
  
  ;loop throught widths to see which is best SNR
  for test_width=1,12 do begin
    test_width_arr=replicate(test_width,n_elements(widtharr))
    
    ;extract new profile
    ;note: using test_width_arr rather than widtharr
    bmep_extraction,j,centerarr,state,order,test_width_arr,bkgndl,bkgndr,$
      printp,fitgaussp,plotp,slidep,pwindowsize,singlep,cosmic_sigma,n_iterate_cosmic,bkgnd_naverage,max_rays_per_col,$
      F,ferr,Fopt,fopterr,img2d_nobkgnd,parr,sky_reslut_arr,sky_residuals,revisevar=state.revisevar,quiet=1
      
    index=where(state.wavel ge leftstatspos and state.wavel le rightstatspos,ct)
    momentresult=moment(f[index],sdev=std_dev)
    momentresultopt=moment(fopt[index],sdev=std_devopt)
    
    ;print info in style of LaTeX
    print,test_width,$ ;width
      ' & ',momentresult[0]/std_dev,$ ; snr box
      ' & ',momentresultopt[0]/std_devopt,$ ; snr optimal
      $;' & ',(momentresultopt[0]/std_devopt)*(momentresult[0]/std_dev),$ ; snr optimal * snr box
      ' & ',momentresult[0],$
      ' & ',momentresultopt[0],$
      ' & ',total(fopt[index]),' \\ \hline' ; mean val
      
  endfor ; test_width
  ;rerun extraction to fix damage.
  bmep_extraction,j,centerarr,state,order,widtharr,bkgndl,bkgndr,$
    printp,fitgaussp,plotp,slidep,pwindowsize,singlep,cosmic_sigma,n_iterate_cosmic,bkgnd_naverage,max_rays_per_col,$
    $;OUTPUTS
    F,ferr,Fopt,fopterr,img2d_nobkgnd,parr,sky_reslut_arr,sky_residuals,revisevar=state.revisevar
  print,'\end{tabular}'
end

;create headers for blindly extracted objects
function bmep_blind_hdr,f,extrainfo1,extrainfo2,yexpect,width,$
    isstar,objnum,min_width,exten=exten,image=image,no_wave_info=no_wave_info
    
  MKHDR, header, f, exten=exten, image=image
  FOR jj=0,n_elements(extrainfo1) -1 do $
    if strnumber(extrainfo2[jj]) then $
    if float(extrainfo2[jj]) eq fix(extrainfo2[jj]) then $
    sxaddpar, Header,extrainfo1[jj],fix(extrainfo2[jj]) else $
    sxaddpar, Header,extrainfo1[jj],float(extrainfo2[jj]) else $
    sxaddpar, Header,extrainfo1[jj],extrainfo2[jj]
  sxaddpar, Header, 'YPOS', yexpect
  sxaddpar, Header, 'WIDTH', width
  sxaddpar, Header, 'OBJNUM', objnum
  sxaddpar, Header, 'ISSTAR', isstar
  sxaddpar, Header, 'WBYHAND', 0
  sxaddpar, Header, 'GAUSRCHI', -1.0
  sxaddpar, Header, 'CBYHAND', 0
  sxaddpar, Header, 'MOM1', -1.0
  sxaddpar, Header, 'MOM2', -1.0
  sxaddpar, Header, 'GAUSSP', 1
  sxaddpar, Header, 'ORDER', -1
  sxaddpar, Header, 'NBINS', 0
  sxaddpar, Header, 'AUTOEX', 1
  sxaddpar, Header, 'BLIND', 1, ' Flag if this extraction was done blindly'
  sxaddpar, Header, 'MINW',min_width
  
  if keyword_set(no_wave_info) then begin
    sxaddpar,Header,'CRVAL1',1
    sxaddpar,Header,'CDELT1',1.0
    sxaddpar,Header,'CRPIX1',1
    sxaddpar,Header,'CTYPE1','LINEAR'
  endif
  
  sxaddpar, Header, 'COMMENT',' Exten 1: Optimal extraction'
  sxaddpar, Header, 'COMMENT',' Exten 2: Optimal extraction error bars'
  sxaddpar, Header, 'COMMENT',' Exten 3: Boxcar extraction'
  sxaddpar, Header, 'COMMENT',' Exten 4: Boxcar extraction error bars'
  sxaddpar, Header, 'COMMENT',' Exten 5: P'
  sxaddpar, Header, 'COMMENT',' This is the blind extraction'
  
  return, header
  
end

;follows horne 86 method of flagging cosmic rays.  it also flags anything
;with a variance of 0 as a bad pixel.
pro bmep_calc_cosmic_rays,i,bkgndl,bkgndr,bottomint,topint,state,totalpixflagged,$
    xarr_big,sky_reslut_arr,p,f,cosmic_sigma,badPixelMask,sky_residuals,n_iterate_cosmic,$
    order,xarr_small,bkgnd_naverage,max_rays_per_col,variance,ferr
  ;variables all the same as in main extraction program.
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  ;all pixels start as good
  badPixelMask=replicate(1.0,n_elements(state.data[0,*]))
  n_flagged=0
  if cosmic_sigma gt 0 then begin
    ;calculate range that you care about.
    if n_elements(bkgndl) ge 1 then begin
      bottom_cosmic_range=fix(min(bkgndl))
      top_cosmic_range=fix(max(bkgndr))
    endif else begin
      bottom_cosmic_range=bottomint
      top_cosmic_range=topint
    endelse
    
    ;check that range is OK
    if bottom_cosmic_range lt 0 then bottom_cosmic_range=0
    if top_cosmic_range gt n_elements(state.data[0,*])-1 then top_cosmic_range=n_elements(state.data[0,*])-1
    if bottom_cosmic_range gt n_elements(state.data[0,*])-1 then bottom_cosmic_range=0
    if top_cosmic_range lt 0  then top_cosmic_range=n_elements(state.data[0,*])-1
    
    ;make names easier...
    bot=bottom_cosmic_range
    top=top_cosmic_range
    
    ;remove zeroes without messing with the removed count
    index=where(variance[i,bot:top] eq 0,ct)
    index=index+bot
    if ct gt 0 then badPixelMask[index]=0.0
    
    ;data - sky - f*p
    residual_arr=(state.data[i,bot:top] - $
      poly(findgen(top-bot+1)+bot,sky_reslut_arr[i,*]) - $
      f[i]*p[bot:top])
    index=where(residual_arr*residual_arr GT $
      cosmic_sigma*cosmic_sigma*variance[i,bot:top],ct)
    index=index+bot
    
    ;actually do the masking if not too many flagged.
    if ct gt 0 and ct lt max_rays_per_col then badPixelMask[index]=0.0 ; mark bad pixels as 0
    nflagged=ct
    
    ;if cosmic rays, iterate things...
    if ct ne 0 and ct ne n_elements(badPixelMask) and ct lt max_rays_per_col then begin
      totalpixflagged=totalpixflagged+ct
      for cosmic_counter=0, n_iterate_cosmic-1 do begin
      
        ;recalculate sky
        newresult=bmep_fit_sky(badPixelMask*state.data[i,*],bkgndl,bkgndr,order,bkgnd_naverage)
        
        ;recalculate boxcar flux
        f[i]=total(badPixelMask*(state.data[i,bottomint:topint] - poly(xarr_small,newresult)))
        ferr[i]=sqrt(total(badPixelMask*(variance[i,bottomint:topint])))
        
        ;recalculate residual array (note new sky, not old)
        residual_arr=badPixelMask[bot:top]*(state.data[i,bot:top] - $
          poly(findgen(top-bot+1)+bot,newresult) - $
          f[i]*p[bot:top])
          
        ;find more bad pixels
        index=where(residual_arr*residual_arr GT cosmic_sigma*cosmic_sigma*variance[i,bot:top],ct) ;flag bad pixels
        
        ;mask pixels
        if ct gt 0 and ct lt max_rays_per_col - nflagged then begin
          badPixelMask[index]=0.0
          totalpixflagged=totalpixflagged+ct
          nflagged=nflagged+ct
        endif
        
        ;bail out of for loop if no moar cosmic rays.
        if ct eq 0 or ct ge max_rays_per_col - nflagged then break
      endfor
      sky_reslut_arr[i,*]=newresult ; update sky array.
    endif
    sky_residuals[i]=total(abs(residual_arr))
  endif;cosmic sigma gt 0
end

;attempt to calculate where the bad skylines are based on variance image
pro bmep_calc_skyline_mask,varimg,percent,rval
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  tempvarimg=varimg
  ;only consider finite values (remove nans)
  index=where(finite(tempvarimg) eq 0,ct)
  tempvarimg[index]=0.0
  
  ;define default return value of 1 for each wavelength
  rval=replicate(1.0,n_elements(varimg[*,0]))
  mask=[]
  
  ;calculate the total noise in one column.
  for i=0,n_elements(varimg[*,0])-1 do mask=[mask,total(varimg[i,*])]
  
  ;calculate bad arrays based on the total noise
  indexbad=where(mask gt bmep_percent_cut(mask,100-percent),maskedct)
  indexgood=where(mask le bmep_percent_cut(mask,100-percent),goodct)
  
  ;0 has unacceptable noise, 1 has good value
  rval[indexbad]=0.0
  rval[indexgood]=1.0
end

;attempt to calculate where the bad skylines are based on noise in extracted spectrum
;so it masks out the noisest lines
pro bmep_calc_skyline_mask_v2,var_opt,percent,rval
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  tempvarimg=var_opt
  index=where(finite(tempvarimg) eq 0,ct)
  tempvarimg[index]=0.0
  
  rval=replicate(1.0,n_elements(var_opt))
  mask=var_opt
  
  indexbad=where(mask gt bmep_percent_cut(mask,100-percent),maskedct)
  indexgood=where(mask le bmep_percent_cut(mask,100-percent),goodct)
  rval[indexbad]=0.0
  rval[indexgood]=1.0
end

;;LEAVE COMMENTED OUT CODE. analytic slitloss in 2d. could be
;;useful if we come back to this method
;function bmep_calc_slitloss_gauss2d,x,y;,x0,y0,xsig,ysig,A
;  restore,'tmp.sav'
;  return,A*exp( -1.0 * ( ( x - x0 )^2 / (2.0*xsig^2) + ( y - y0 )^2 / (2.0*ysig^2) ) )
;end
;
;;LEAVE COMMENTED OUT CODE. analytic slitloss in 2d. could be
;;useful if we come back to this method
;function bmep_calc_slitloss_ylimits,y
;  restore,'tmp.sav'
;  return,[0.0,ysig*n_sigma_width*2.0]
;end

;input: gauss_sigma. sigma of the gaussian fit to the y-profile.
function bmep_calc_slitloss,gauss_sigma
  pixscale=0.1799 ;pixelscale in arcsec/pixel
  slitwidth=0.7 ; slitwidth in arcseconds
  return,erf((slitwidth/pixscale/2)/(gauss_sigma)/sqrt(2.0))
;
;
;
;;LEAVE COMMENTED OUT CODE. analytic slitloss in 2d. could be
;;useful if we come back to this method
;  n_sigma_width=6.0 ; number of sigma to make
;    ;so for the "6" the total integrated area
;    ;is plus minus 6 sigma from the center.
;  xsig=gauss_sigma;/(2.0*SQRT(2.0*ALOG(2.0)))
;  ysig=xsig
;  x0=xsig*n_sigma_width ;put x0 in the center
;  y0=ysig*n_sigma_width ;put y0 in the center
;  A=1.0 ; arbitrary constant.
;  save,x0,y0,xsig,ysig,A,n_sigma_width,filename='tmp.sav'
;  all_flux=int_2d('bmep_calc_slitloss_gauss2d',$
;    [0.0,xsig*n_sigma_width*2.0],'bmep_calc_slitloss_ylimits',96,/double)
;
;  ;integrate a 0.7 slit across gaussian
;  midpoint=xsig*n_sigma_width
;  seen_flux=int_2d('bmep_calc_slitloss_gauss2d',$
;    [midpoint-((slitwidth/pixscale)/2.0),midpoint+((slitwidth/pixscale)/2.0)],$
;    'bmep_calc_slitloss_ylimits',96,/double)
;  print,all_flux
;  print,seen_flux
;  spawn,'rm tmp.sav'
;return,seen_flux/all_flux
end


;check if a directory exists.
;stolen from somewhere....
function bmep_dir_exist, dir
  CD, CUR=cur
  CATCH, error_status
  if (error_status NE 0) then begin
    return, 0
  endif
  CD, dir
  CD, cur
  return, 1
end


; Anyone can run BMEP on any normal fits image by calling this function
; This program calls the actual program, bmep and bmep_lris are instrument
; specific wrappers.  BMEP is for mosfire, bmep_lris is for LRIS data.
;
; Inputs
;big_img - image to be displayed. can be different than science image if you
;            wanted to show the variance image as well.  The X SIZE of this image
;            should be the same as the science image
;
;sciimg - science image in number of electrons.
;
;var_img - variance of the above image.  If assuming poisson noise, this input would
;          be exactly the same as the science image
;
;highval - Sometimes there are cosmic rays that mess up the scale of the images and
;          make it impossible to see the image properly.  If there is one spike at,
;          say, a value of 60000 but the majority of your image is between 3000 and
;          5000, you would want to set this value to something like 5500. This is because
;          there is basically no data between 5000 and 60000 but this eats up all
;          the dynamic range of the image display. to exclude the top 10 percent of
;          the pixels (set them to white...) make this bmep_percent_cut(sciimg,90).
;          To include the full range of the image set this to max(sciimg).
;
;lowval - Same thing as highval but this sets the lower limit.  To exclude
;
;slitname  - name of object in slit for saving purposes.
;
;filtername - name of filter for slit for saving purposes.
;
;wavel - wavelength array of the data.  Should have same number of elements as x direction
;        of science image.
;
;savepath - output will be saved in savepath+slitname+'.sav'
;           to save in current directory (not advisable) use set savepathe to ''
;           if multiple objects are in slit and they are extracted, the path will be
;           savepath+slitname+suffix+'.sav' where "suffix" is "--N" where N is the
;           object number.
;
;revisevar - set keyword to revise variance estimate. no default, please set this. (essential keyword)
;
;extrainfo1 - fits keywords to be added to fits header (optional)
;extrainfo2 - values of the fits keywords to be added (optional)
;extrainfo2 - comments of the fits keywords to be added (optional)
;savetext - flag as to weather or not to save a text file. (optional)
;
;
;NOTES: this will only work on idl 8.1 or higher.
;       The majority of the program is housed in the bmep_keyboardhandler() function.
;


pro bmep_display_image,big_img,sciimg,var_img,highval,lowval,slitname,filtername,wavel,savepath,$
    revisevar=revisevar,extrainfo1=extrainfo1,extrainfo2=extrainfo2,extrainfo3=extrainfo3,savetext=savetext
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  if ~keyword_set(extrainfo1) then extrainfo1=['0']
  if ~keyword_set(extrainfo2) then extrainfo2=['0']
  if ~keyword_set(extrainfo3) then extrainfo3=['0']
  if ~keyword_set(savetext) then savetext=0
  !except=2
  
  ;check version number
  IF (Float(!Version.Release) lt 8.1) THEN message,'idl ver must be 8.1 or greater'
  
  index=where(var_img eq 0,ct)
  ;print,'there are ',ct,' pixels with 0 variance '
  
  index=where(finite(var_img) eq 0,ct)
  ;print,'there are ',ct,' pixels with nan value for  variance '
  sciimg[index]=0.0
  var_img[index]=0.0
  
  index=where(var_img lt 0,ct)
  var_img[INDEX]=ABS(var_img[INDEX])
  ;print,'there are ',ct,' pixels with NEGATIVE value for  variance '
  
  im1=image(big_img,dimensions=[n_elements(big_img[*,0]),n_elements(big_img[0,*])],$
    window_title=slitname+'-'+filtername,max_value=highval,min_value=lowval,margin=[0,0,0,0],$
    location=[0,0],/xstyle,/ystyle);title=slitname,
    
  im1.window.UVALUE={x0:0, y0:0, buttonDown:0L,$
    ydata:replicate(0.D,n_elements(sciimg[0,*])),ydataerr:replicate(0.D,n_elements(sciimg[0,*])),data:sciimg,$
    slitname:slitname+'_'+filtername,wavel:wavel,savepath:savepath,var_img:var_img,revisevar:revisevar,$
    raw_slitname:slitname,extrainfo1:extrainfo1,extrainfo2:extrainfo2,extrainfo3:extrainfo3,savetext:savetext,cont_mode:0,$
    stats_mode:0,n_bins:0,l_bins:[0,0,0,0,0,0,0,0,0,0,0],r_bins:[0,0,0,0,0,0,0,0,0,0,0]}
  im1.window.MOUSE_DOWN_HANDLER='bmep_MouseDown'
  im1.window.MOUSE_UP_HANDLER='bmep_MouseUp'
  im1.window.Keyboard_Handler='bmep_KeyboardHandler'
  
end


;extract the spectra!! algorithm from horne 1986
pro bmep_extraction_simple,sciimg,var_img,ydata,ypos,width,f,ferr,fopt,fopterr,p

  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
    
  cosmic_sigma=4.5
  max_rays_per_col=10
  f=[]
  ferr=[]
  fopt=[]
  fopterr=[]
  
  bottomint=fix(ypos-width)
  if bottomint lt 0 then bottomint=0
  topint=fix(ypos+width)
  if topint gt n_elements(sciimg[0,*])-1 then topint=n_elements(sciimg[0,*])-1
  
  xarr_small=findgen(topint-bottomint+1)+bottomint
  xarr_big=findgen(n_elements(ydata))
  
  ;double check that the extraction ranges are actually the same.
  if min(xarr_small) ne bottomint then print,'WARNING, THE EXTRACTION RANGES ARE DIFFERENT.'
  variance=var_img
  for i=0,n_elements(sciimg[*,0])-1 do begin
    f=[f,total(sciimg[i,bottomint:topint])]
    ferr=[ferr,sqrt(total(variance[i,bottomint:topint]))]
  endfor
  
  ;===OPTIMAL===================================================
  p=ydata
  index=where(xarr_big lt bottomint or xarr_big gt topint,/null)
  p[index]=0.0
  index=where(p lt 0,/null)
  p[index]=0.0
  if total(p) ne 0 then p=p/total(p)
  gresult=MPFITPEAK(double(xarr_big),double(p),$
    coeff,nterms=3,/gaussian)
  p=gresult
  index=where(xarr_big lt bottomint or xarr_big gt topint,ct)
  if ct ne 0 then p[index]=0.0
  index=where(p lt 0,ct)
  if ct gt 0 then p[index]=0.0
  if total(p) ne 0 then p=p/total(p)
  
  
  for i=0,n_elements(sciimg[*,0])-1 do begin
  
    badPixelMask=replicate(1.0,n_elements(p))
    residual_arr=(sciimg[i,xarr_small] - f[i]*p[xarr_small])
    index=where(residual_arr*residual_arr GT $
      cosmic_sigma*cosmic_sigma*variance[i,xarr_small],ct)
    if ct lt max_rays_per_col and ct gt 0 then badPixelMask[xarr_small[index]]=0.0
    
    var_opt=variance[i,*]
    
    index=where(var_opt[xarr_small] eq 0.0 or var_opt[xarr_small] eq 99.00 or var_opt[xarr_small] eq 99.0*99.0,count)
    if count gt 0 then begin
      badPixelMask[xarr_small[index]]=0.0
      var_opt[xarr_small[index]]=1.0
    endif ;else begin
    
    ;calc optimal flux as in horne
    numerator=total(badPixelMask[xarr_small]*p[xarr_small]*$
      (sciimg[i,xarr_small])/var_opt[xarr_small])
    denominator=total(badPixelMask[xarr_small]*p[xarr_small]*p[xarr_small]/var_opt[xarr_small])
    if denominator ne 0 then flux_opt=numerator/denominator else flux_opt=0.0
    
    ;calc optimal flux error as in horne.
    numerator=total(badPixelMask[xarr_small]*p[xarr_small]) ;
    denominator=total(badPixelMask[xarr_small]*p[xarr_small]*p[xarr_small]/var_opt[xarr_small])
    if denominator ne 0 then err_opt=sqrt(abs(numerator/denominator)) else err_opt=0.0
    
    ;undo the damage from setting pixels with zero variance to one
    if count gt 0 then var_opt[xarr_small[index]]=0.0
    
    fopt=[fopt,flux_opt] ;!!!!!!!!!!!!
    fopterr=[fopterr,err_opt]
  endfor ; cols of data! - end of the optimal section of extraction.
  
end ; end of extraction


;extract the spectra!! algorithm from horne 1986
pro bmep_extraction,j,centerarr,state,order,widtharr,bkgndl,bkgndr,$
    printp,fitgaussp,plotp,slidep,pwindowsize,singlep,cosmic_sigma,n_iterate_cosmic,bkgnd_naverage,max_rays_per_col,$
    $;OUTPUTS
    F,ferr,Fopt,fopterr,img2d_nobkgnd,parr,sky_reslut_arr,sky_residuals,$
    revisevar=revisevar,quiet=quiet
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
    
  if not keyword_set(quiet) then quiet = 0
  
  ;print,'doing center number ',j
  f=[]
  ferr=[]
  fopt=[]
  fopterr=[]
  sky_reslut_arr=findgen(n_elements(state.data[*,0]),order+1)
  
  ;lower extraction limit
  bottomint=fix(centerarr[j]-widtharr[j])
  if bottomint lt 0 then bottomint=0
  
  ;upper extraction limit
  topint=fix(centerarr[j]+widtharr[j])
  if topint gt n_elements(state.data[0,*])-1 then topint=n_elements(state.data[0,*])-1
  
  ;calculate array to extract
  xarr_small=findgen(topint-bottomint+1)+bottomint
  
  ;calculate array for whole column
  xarr_big=findgen(n_elements(state.ydata))
  
  ;double check that the extraction ranges are actually the same.
  if min(xarr_small) ne bottomint then print,'WARNING, THE EXTRACTION RANGES ARE DIFFERENT.'
  if max(xarr_small) ne topint    then print,'WARNING, THE EXTRACTION RANGES ARE DIFFERENT.'
  ;PRINT,'Extracting from ',bottomint,' to ',topint,'    ',n_elements(xarr_small),' elements'
  
  ;STEP2 (steps from horne 1986 extraction paper... step 1 was flatfielding)
  variance=state.var_img
  skysub_image=replicate(0.0,n_elements(state.data[*,0]),n_elements(state.data[0,*]))
  for i=0,n_elements(state.data[*,0])-1 do begin
  
    ;STEP3  (calculate sky)
    result=bmep_fit_sky(state.data[i,*],bkgndl,bkgndr,order,bkgnd_naverage)
    SKY=poly(xarr_small,result)
    sky_reslut_arr[i,*]=result
    
    ;step4 (extract spectrum and error)
    f=[f,total(state.data[i,bottomint:topint] - sky)]
    ferr=[ferr,sqrt(total(variance[i,bottomint:topint]))]
    
  endfor ; looping through i, the number of columns...
  
  ;===OPTIMAL===================================================
  ;optimal extraction based on horne 1989
  ;step 5 (calculate p)
  
  ;figure out p from cross section and use only that.
  if singlep then begin
    if slidep then message,'insanity occured. only single p or slide p allowed'
    
    FORWARD_FUNCTION bmep_find_p
    ;calculate p
    p=bmep_find_p(state,bkgndl,bkgndr,order,bottomint,topint,printp,fitgaussp,plotp,xarr_big,bkgnd_naverage)
    
    ;plot p on the plot
    plotshift=poly(centerarr[j],bmep_fit_sky(state.ydata,bkgndl,bkgndr,order,bkgnd_naverage))
    pscale=(state.ydata[centerarr[j]]-plotshift)/max(p)
    oplot,(p*pscale + plotshift),color=255,linestyle=3
  endif
  
  ;set up vars
  parr=replicate(0.0,n_elements(state.data[*,0]),n_elements(state.data[0,*]))
  img2d_nobkgnd=replicate(0.0,n_elements(state.data[*,0]),n_elements(state.data[0,*]))
  sky_residuals=replicate(0.0,n_elements(state.data[*,0]))
  totalpixflagged=0
  
  ;loop through columns
  for i=0,n_elements(state.data[*,0])-1 do begin
  
    ;if spectrum is UNDERSAMPLED in the spatial direction,
    ;p (the fraction of light in each pixel) could change as
    ;a function of wavelength.  If this is true, p should be
    ;reevaluated as a function of position. The way I do this
    ;pretty much only works for anything bright.
    ;(step 5) Alternate version of step 5...
    if slidep then begin
      if singlep then message,'insanity occured. only single p or slide p allowed'
      p=bmep_find_p_slide(state,pwindowsize,i,f,xarr_small)
    endif
    
    ;step 7 (mask cosmic rays)
    ;calculate range over which we care about cosmic rays
    if keyword_set(revisevar) then bmep_calc_cosmic_rays,i,bkgndl,bkgndr,bottomint,topint,state,totalpixflagged,$
      xarr_big,sky_reslut_arr,p,f,cosmic_sigma,badPixelMask,sky_residuals,n_iterate_cosmic,$
      order,xarr_small,bkgnd_naverage,max_rays_per_col,variance,ferr $
    else badPixelMask=replicate(1.0,n_elements(xarr_big))
    
    ;swapped steps 6 and 7 to use better flux estimation after removing cosmic spikes.
    ;step 6 (Revise errors) (dont revise for MOSFIRE)
    if keyword_set(revisevar) then var_opt=poly(xarr_big,sky_reslut_arr[i,*]) + f[i]*p $
    else var_opt=variance[i,*]
    
    ;step 8
    ;try to not divide by 0.
    index=where(var_opt[xarr_small] eq 0.0 or var_opt[xarr_small] eq 99.00,count)
    if count gt 0 then begin
      ;          print,i,count,' number of var_opt bade '
      ;          print,index
      ;places with 0 variance are considered bad
      badPixelMask[xarr_small[index]]=0.0
      ;though they are still in the denominator, so set
      ;their value to 1.0 temporarily.
      ;these points are masked, so their value doesn't matter
      var_opt[xarr_small[index]]=1.0
    endif ;else begin
    
    ;calc optimal flux as in horne
    numerator=total(badPixelMask[xarr_small]*p[xarr_small]*$
      (state.data[i,xarr_small] - poly(xarr_small,sky_reslut_arr[i,*]))/var_opt[xarr_small])
    denominator=total(badPixelMask[xarr_small]*p[xarr_small]*p[xarr_small]/var_opt[xarr_small])
    if denominator ne 0 then flux_opt=numerator/denominator else flux_opt=0.0
    
    ;calc optimal flux error as in horne.
    numerator=total(badPixelMask[xarr_small]*p[xarr_small]) ;
    denominator=total(badPixelMask[xarr_small]*p[xarr_small]*p[xarr_small]/var_opt[xarr_small])
    if denominator ne 0 then err_opt=sqrt(abs(numerator/denominator)) else err_opt=0.0
    
    ;undo the damage from setting pixels with zero variance to one
    if count gt 0 then var_opt[xarr_small[index]]=0.0
    
    fopt=[fopt,flux_opt] ;!!!!!!!!!!!!
    fopterr=[fopterr,err_opt]
    
    img2d_nobkgnd[i,*]=state.data[i,*] - poly(findgen(n_elements(state.data[i,*])),sky_reslut_arr[i,*])
    parr[i,*]=p
  ;        if i gt 405 and flux_opt eq 0 then stop
  endfor ; cols of data! - end of the optimal section of extraction.
  if not quiet then print,'eliminated cosmic rays: ',totalpixflagged
end ; end of extraction




;
;fit the sky background to a polynomial.
;INPUTS:
;  yvals- y values of data cross-section
;  bkgndl - array of left background values
;  bkgndr - array of right background values
;  order - order of poly to be fit
;  bkgnd_naverage - how many pixels to average in
;           doing the background. just like IRAF b_naverage
;  overplotflag- it will overplot the background on current
;      plotting device (1 to plot, 0 to not)
;
;OUTPUT:
;  return value - result of poly_fit() on the background.
;      Will return zeroes if something happened (like not
;      enough points...)
;
;DESCRIPTION:
;  function to fit the sky using a poly fit.
;  0 is considered to be a bad pixel and isn't used in the fit
;  although 0 *could* be a prefectly
;  legitimate value for a pixel,
;  it is EXTREMELY unlikely that a pixel
;  will ever be EXACTLY the value of 0
;
function bmep_fit_sky,yvals,bkgndl,bkgndr,order,bkgnd_naverage,overplotflag=overplotflag
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  ;Case of no background return line at y=0
  if n_elements(bkgndl) eq 0 then return,replicate(0.0,order+1)
  
  xarr_big=findgen(n_elements(yvals))
  ;Find bad pixels
  index=where(yvals ne 0,ct)
  
  ;Case of no good pixels
  if ct eq 0 then return,replicate(0.0,order+1)
  
  ;Mask bad pixels
  yvals=yvals[index]
  xarr_big=xarr_big[index]
  
  ;determine the x and y values of the data to fit
  ;binning if appropraite.
  xfit=[] ; x data to fit
  yfit=[] ; y data to fit
  for k=0,n_elements(bkgndl)-1 do begin
    index=where(xarr_big ge bkgndl[k] and xarr_big le bkgndr[k],ct)
    if ct gt 0 then begin
      if abs(bkgnd_naverage) gt 1 then begin
        if n_elements(xarr_big[index]) ge abs(bkgnd_naverage) then begin
          summed_yvals=[]
          summed_xvals=[]
          i=0
          while i lt n_elements(index) - abs(bkgnd_naverage) +1 do begin
            if bkgnd_naverage lt 0 then summed_yvals=[summed_yvals,median(yvals[index[i:i+abs(bkgnd_naverage)-1]])]
            if bkgnd_naverage gt 0 then summed_yvals=[summed_yvals,avg(yvals[index[i:i+abs(bkgnd_naverage)-1]])]
            summed_xvals=[summed_xvals,(float(xarr_big[index[i]]) + float(xarr_big[index[i+abs(bkgnd_naverage)-1]]))/2.0]
            i=i+abs(bkgnd_naverage) ;+1
          endwhile
          xfit=[xfit,summed_xvals]
          yfit=[yfit,summed_yvals]
        endif else begin ; too few points to avg.
          xfit=[xfit,xarr_big[index]]
          yfit=[yfit,yvals[index]]
        endelse
      endif else begin ;if no averaging/median of bkgnd points
        xfit=[xfit,xarr_big[index]]
        yfit=[yfit,yvals[index]]
      endelse
    endif;ct gt 0
  endfor ; each background spot
  
  ;Case of no good data to use for bkgnd
  if n_elements(xfit) eq 0 then return,replicate(0.0,order+1)
  
  ;Remove duplicates
  index=rem_dup(xfit)
  xfit=xfit[index]
  yfit=yfit[index]
  
  if N_elements(xfit) le order then return,replicate(0.0,order+1) else $
    result=poly_fit(xfit,yfit,order,yfit=output,status=poly_status,MEASURE_ERRORS=sqrt(yfit),sigma=resulterror)
    
  ;Check for an error in the fit
  if poly_status eq 1 then begin
    print,'no good pixels in the sky?!?'
    return,replicate(0.0,order+1)
  endif
  
  if keyword_set(overplotflag) then oplot,xfit,yfit,psym=6,thick=2
  return,result
end


;calculate p as a function of position using a small window
;of size 2*windowsize
function bmep_find_p_slide,state,pwindowsize,i,f,xarr_small,sky_reslut_arr
  p=replicate(0.0,n_elements(state.data[0,*]))
  
  ;create binned data to measure p
  if i - pwindowsize lt 0 then left=0 else left = i - pwindowsize
  if i + pwindowsize gt n_elements(state.data[*,0])-1 then $
    right=n_elements(state.data[*,0])-1 else right = i + pwindowsize
    
  ;calc p
  for counter=left,right do $
    if f[counter] ne 0 then $
    p[xarr_small]=p[xarr_small]+(state.data[counter,xarr_small]-poly(xarr_small,sky_reslut_arr[counter,*]))/f[counter]
    
  ;normalization and negative elemination
  index=where(p lt 0,ct)
  if ct gt 0 then p[index]=0.0
  if total(p) ne 0 then p=p/total(p)
  return,p
end

;;calculate p as a function of position using a small window
;;of size 2*windowsize
;function bmep_find_p_slide_v2,state,xarr_small,sky_reslut_arr,centerarr,widtharr,f,j
;
;  ;lower extraction limit
;  bottomint=fix(centerarr[j]-widtharr[j])
;  if bottomint lt 0 then bottomint=0
;
;  ;upper extraction limit
;  topint=fix(centerarr[j]+widtharr[j])
;  if topint gt n_elements(state.data[0,*])-1 then topint=n_elements(state.data[0,*])-1
;  parr=[]
;  skysub_img=state.data
;
;   x_arr_big=findgen(n_elements(state.data[0,*]))
;
;   ;subtract the sky and normalize by flux
;   for i=0,n_elements(state.data[*,0])-1 do begin
;     skysub_img[i,*]= (state.data[i,*]-poly(x_arr_big,sky_reslut_arr[i,*]))
;     deno=total((state.data[i,*]-poly(x_arr_big,sky_reslut_arr[i,*])))
;     if deno ne 0 then skysub_img[i,*]=skysub_img[i,*]/deno
;     index=where(skysub_img[i,*] gt 1 or skysub_img[i,*] lt 0,ct)
;     if ct ne 0 then skysub_img[i,index]=0
;     endfor
;   ;make outsides zero.
;   !p.multi=[0,1,1]
;   plot,skysub_img[*,bottomint+(topint-bottomint)/2],nsum=5
;   oplot,skysub_img[*,bottomint+(topint-bottomint)/2+2],nsum=5,color=254
;   im1=image(skysub_img,dimensions=[n_elements(skysub_img[*,0]),n_elements(skysub_img[0,*])])
;   stop
;   for i=0,n_elements(state.data[*,0])-1 do begin
;    p=replicate(0.0,n_elements(state.data[0,*]))
;
;    endfor
;  return,parr
;end

;calculate spatial profile of the object's light.
function bmep_find_p,state,bkgndl,bkgndr,order,bottomint,topint,printp,fitgaussp,plotp,xarr_big,bkgnd_naverage

  ;calculate sky
  result=bmep_fit_sky(state.ydata,bkgndl,bkgndr,order,bkgnd_naverage)
  
  ;calculate p
  p=state.ydata-poly(xarr_big,result)
  ;  p=smooth(p,3,/edge_truncate) ; DON'T SMOOTH THAT SHIT!!
  
  ; outside of width, p is zero
  index=where(xarr_big lt bottomint or xarr_big gt topint,ct)
  if ct ne 0 then p[index]=0.0
  
  ;enforce positivity
  index=where(p lt 0,ct)
  if ct gt 0 then p[index]=0.0
  
  ;normalize
  if total(p) ne 0 then p=p/total(p)
  
  ;print if desired
  if printp then print,'p'
  if printp then print,p
  if printp then print
  
  ;fit p profile to a gaussian if desired
  if fitgaussp then begin
    ;gresult=gaussfit(double(xarr_big),double(p),nterms=3)
    gresult=MPFITPEAK(double(xarr_big),double(p),$
      coeff,nterms=3,/gaussian)
    p=gresult
    
    ;set p to 0 outside the width of extraction.
    index=where(xarr_big lt bottomint or xarr_big gt topint,ct)
    if ct ne 0 then p[index]=0.0
    
    ;normalize p
    index=where(p lt 0,ct)
    if ct gt 0 then p[index]=0.0
    if total(p) ne 0 then p=p/total(p)
    
    
    ;print if desired
    if printp then print,'pgauss'
    if printp then print,p
    if printp then print
  endif
  
  ;plot if desired
  if plotp then zzzzzz=plot(p/max(p))
  
  ;return
  return,p
end




;input a filename and extract a slit if in form of
;  maskname_fliter_slitname_ivar.fits.gz
;  or
;  slitname_fliter_eps.fits.gz
;  ex:
;  cosmos_z3_01_H_1823_ivar.fits.gz
;
;and for lris
;  cutout_slitname.fits
;and program will extract only the slitname
function bmep_get_slitname,filename,maskname,ivar=ivar,eps=eps,lris=lris,idlreduc=idlreduc,gzending=gzending
  postlength=0
  if keyword_set(lris) then begin
    start=strlen('cutout_')
    postlength=strlen(filename)-start-5 ;  5 for the .fits
    return,strcompress((strmid(filename,start,postlength)),/remove_all)
  endif ;lris data...
  if keyword_set(ivar) then begin
    postlength=strlen(filename)-strlen('_ivar.fits.gz') - (strlen(maskname)+3)
    if ~keyword_set(gzending) then postlength=strlen(filename)-strlen('_ivar.fits') - (strlen(maskname)+3)
    return,strcompress((strmid(filename,strlen(maskname)+3,postlength)),/remove_all)
  endif
  if keyword_set(eps) then begin
    postlength=strlen(filename)-strlen('_eps.fits.gz') - (strlen(maskname)+3)
    if ~keyword_set(gzending) then postlength=strlen(filename)-strlen('_eps.fits') - (strlen(maskname)+3)
    return,strcompress((strmid(filename,strlen(maskname)+3,postlength)),/remove_all)
  endif
  if keyword_set(idlreduc) then begin
    postlength=strlen(filename)-strlen('.2d.fits') -strlen('co2_03.')
    return,strcompress((strmid(filename,strlen('co2_03.'),postlength)),/remove_all)
  endif
  return,'NO_SLITNAME_FOUND'
end



function bmep_guess_redshift,state,wavel,fopt,linename=linename,linecol=linecol
  ;get mask name
  index=where(state.extrainfo1 eq 'MSKNM',ct)
  if ct eq 0 then begin
    linename='none'
    linecol=-999
    return,-999.9
  endif
  maskname=state.extrainfo2[index[0]]
  z=strmid(maskname,2,1)
  
  ;get filtername
  index=where(state.extrainfo1 eq 'FILTNM')
  filtername=strcompress(state.extrainfo2[index],/remove_all)
  
  ;define line names
  linenames=['[OII]', 'Hb',  '[OIII]','Ha']
  linewavels=[3728.5,4862.7,5008.2,6564.6]
  
  best_ind=0
  ;get bins
  NBINS=state.n_bins
  if NBINS gt 1 then begin
    tot_fluxes=[]
    for i=0,NBINS-1 do begin
      tot_fluxes=[tot_fluxes,total(fopt[state.l_bins[i]:state.r_bins[i]])]
      if abs(state.l_bins[i]-state.r_bins[i]) gt 100 then tot_fluxes[-1]=0.0
    endfor
    best_ind=where(tot_fluxes eq max(tot_fluxes),ct)
    best_ind=best_ind[0]
  endif
  
  LBIN1=state.l_bins[best_ind]
  RBIN1=state.r_bins[best_ind]
  
  
  temp=max(fopt[LBIN1:RBIN1],wsub)
  Wavel_guess=state.wavel[wsub]
  
  ;check for case of only wide selections
  if n_elements(fopt[LBIN1:RBIN1]) gt 100 then begin
    linecol=-1
    linename='???'
    redshift=-999.9
    print,'WARNING: ONLY WIDE SELECTIONS. NO Z DETERMINABLE'
    print,'filtername ',filtername
    print,'linename ','???'
    print,'lineguess ',lineguess
    print,'redshift ',redshift
    print,'Wavel_guess ',0
    print,'linecol ',linecol
    return,redshift
  endif
  
  lineguess=0
  CASE z OF
    '1': begin
      if filtername eq 'Y' then begin
        linename='[OII]'
        lineguess=linewavels[where(linenames eq '[OII]')]
      endif else if filtername eq 'J' then begin
        linename='OIII5008'
        lineguess=linewavels[where(linenames eq '[OIII]')]
      endif else if filtername eq 'H' then begin
        linename='HA'
        lineguess=linewavels[where(linenames eq 'Ha')]
      endif else if filtername eq 'K' then begin
        linename='??'
        lineguess=0.0
      endif
    end
    '2': begin
      if filtername eq 'Y' then begin
        linename='??'
        lineguess=0.0
      endif else if filtername eq 'J' then begin
        linename='OII'
        lineguess=linewavels[where(linenames eq '[OII]')]
      endif else if filtername eq 'H' then begin
        linename='OIII5008'
        lineguess=linewavels[where(linenames eq '[OIII]')]
      endif else if filtername eq 'K' then begin
        linename='HA'
        lineguess=linewavels[where(linenames eq 'Ha')]
      endif
    end
    '3': begin
      if filtername eq 'Y' then begin
        linename='??'
        lineguess=0.0
      endif else if filtername eq 'J' then begin
        linename='??'
        lineguess=0.0
      endif else if filtername eq 'H' then begin
        linename='OII'
        lineguess=linewavels[where(linenames eq '[OII]')]
      endif else if filtername eq 'K' then begin
        linename='OIII5008'
        lineguess=linewavels[where(linenames eq '[OIII]')]
      endif
    end
  endcase
  lineguess=lineguess[0]
  
  if lineguess ne 0 then redshift=(Wavel_guess-lineguess)/lineguess else redshift=0.0
  linecol=round(avg([LBIN1,RBIN1]))
  print,'z',z
  print,'linename ',linename
  print,'lineguess ',lineguess
  print,'redshift ',redshift
  print,'wavel guess ',Wavel_guess
  print,'bin number ',best_ind
  
  return,redshift
  
end



FUNCTION bmep_KeyboardHandler, oWin, $
    IsASCII, Character, KeyValue, X, Y, Press, Release, KeyMods
    
  state = oWin.uvalue
  IF Release then return, 1
  
  CASE string(character) OF
    'b': begin
      print,'doing the background, press q to exit'
      order=1 ; order of bkgnd fit
      default_width=4 ; width of extraction
      
      ;cosmic ray settings
      n_iterate_cosmic=3 ; max possible iterations of cosmic ray rejections, like iraf (3)
      cosmic_sigma=4.5 ; default sigma threshold for CR rejection (4.5)
      max_rays_per_col=2 ; maximum number of cosmic rays per column rejected (2/3 for mosfire,
      ; 10 or so for LRIS or stuff with sky)
      
      ;sky settings
      bkgnd_naverage=-3 ; like b_naverage in iraf. (-3)
      ; instead of using each pixel in the background, use an
      ; average or median (negative value for median)
      skymaskpercent=00.0 ; do not show extracted flux that is in the worst N percent
      ; where N is the value of this parameter.  In the extracted
      ;spectrum, the masked pixels are shown as red points and are
      ;not connected to the other points.  It really makes the plot
      ;look nice. (0.00 or ~15 to 25)
      showskymask=0 ;display 2d image showing where the sky has been masked.
      autosky=0 ; if you use "m" to get the center, then automatically select the sky (0 for mosfire, 1 for LRIS)
      autoskybuffer=3 ; distance to the left and right of edge of extraction to start sky (3)
      autoskywidth=12 ; width of each side of the sky.(12)
      mwidth=8 ; width used to get FWHM of profile.(8)
      
      ;p settings
      
      ;choose one, slide p or single p (recommended settings in ())
      singlep=1 ;flag -calc p once based on selected profile (1)
      slidep=0 ;flag  -calc p as a function of position using a window of pixels around every pixel (0)
      pwindowsize=n_elements(state.data[*,0])*0.10/2 ; window size (in num pixels) (1/10th of your npixels)
      fitgaussp=1 ;fit p to a gaussian? (only for calc p once) (1)
      printp=0 ;flag (0)
      plotp=0 ;flag  (0)
      
      
      ;default values, not settings (no change plz)
      objnum=-1
      xr=[min(state.wavel),max(state.wavel)]
      key=' '
      ;each index is a different extraction for the following:
      bkgndl=[]
      bkgndr=[]
      centerarr=[]
      gausschiarr=[]
      wbyhand=[]
      cbyhand=[]
      widtharr=[]
      mom1=[]
      mom2=[]
      slitloss=[]
      gwidth=[]
      gcenter=[]
      gamp=[]
      glinear=[]
      gwidth_err=[]
      gcenter_err=[]
      gamp_err=[]
      glinear_err=[]
      
      
      ;other settings.
      extraction=0 ;flag
      leftstatspos=-1.0
      rightstatspos=-1.0
      saveflag=0  ;flag
      overploterr=0 ;flag
      findautowidth=1 ;flag
      autoextractflag=1 ;flag
      
      ;normalize cause of errors.
      state.ydataerr=sqrt(state.ydataerr)/max(state.ydata) ; convert error to stddev
      state.ydata=state.ydata/max(state.ydata)
      
      ;;uncomment for LRIS settings
      ;skymaskpercent=0.0
      while key ne 'q' and key ne "Q" do begin
        statshift=0.0
        ;plot the things
        if extraction then !p.multi=[0,1,n_elements(centerarr)*2+1] else !p.multi=[0,1,1]
        plot,state.ydata,/ynozero,title=state.slitname
        for i=0,n_elements(bkgndl)-1 do oplot_vert,bkgndl[i],minmax(state.ydata),2
        for i=0,n_elements(bkgndr)-1 do oplot_vert,bkgndr[i],minmax(state.ydata),2
        
        ;if center of profile was found, plot it
        if n_elements(centerarr) ne 0 then $
          for i=0,n_elements(centerarr)-1 do begin
          oplot_vert,centerarr[i],minmax(state.ydata),3
          oplot_vert,centerarr[i]-widtharr[i],minmax(state.ydata),3
          oplot_vert,centerarr[i]+widtharr[i],minmax(state.ydata),3
        ;            print,'plot vert min: ',centerarr[i]-widtharr[i]
        ;            print,'plot vert max: ',centerarr[i]+widtharr[i]
        endfor
        ;find the fitted background
        result=bmep_fit_sky(state.ydata,bkgndl,bkgndr,order,bkgnd_naverage,overplotflag=1)
        ;plot the fitted background
        oplot,poly(findgen(n_elements(state.ydata)),result),linestyle=3,color=255
        
        ;do the extraction
        if extraction then begin
          for j=0, n_elements(centerarr)-1 do begin
            ;extraction!!
            bmep_extraction,j,centerarr,state,order,widtharr,bkgndl,bkgndr,$
              printp,fitgaussp,plotp,slidep,pwindowsize,singlep,cosmic_sigma,n_iterate_cosmic,bkgnd_naverage,max_rays_per_col,$
              $;OUTPUTS
              F,ferr,Fopt,fopterr,img2d_nobkgnd,parr,sky_reslut_arr,sky_residuals,revisevar=state.revisevar
              
              
            ;                if state.savetext eq 0 then begin
            ;                ;plot background subtracted image.
            ;                if n_elements(im) ne 0 then oplotval=im else oplotval=0
            ;;                im=image(img2d_nobkgnd,dimensions=[n_elements(img2d_nobkgnd[*,0]),$
            ;;                n_elements(img2d_nobkgnd[0,*])],overplot=oplotval,$
            ;;                  min_value=bmep_percent_cut(abs(img2d_nobkgnd),5.0),$
            ;;                  max_value=bmep_percent_cut(abs(img2d_nobkgnd),95.0),location=[0,0],margin=[0,0,0,0])
            ;                endif ;save txt? (view bkgnd subbed img)
              
            ;do not plot bad pixels
            bmep_calc_skyline_mask_v2,ferr,skymaskpercent,skymask
            
            yr=[bmep_percent_cut(fopt,2),bmep_percent_cut(fopt,98)]
            index=where(fopt ge min(xr) and fopt le max(xr),ct)
            if ct gt 5 and ct lt 300 then yr=minmax(fopt[index])
            
            plot,state.wavel,f*skymask,xrange=xr,$
              yrange=yr,/nodata
            oplot,state.wavel,f,color=254,psym=3
            oplot,state.wavel,f*skymask
            if overploterr then oploterr,state.wavel,f,ferr,3
            NBINS=state.n_bins
            for k=0,NBINS-1 do begin
              oplot,[state.wavel[state.l_bins[k]],state.wavel[state.l_bins[k]]],minmax(f),color=255
              oplot,[state.wavel[state.r_bins[k]],state.wavel[state.r_bins[k]]],minmax(f),color=255
              xyouts,avg([state.wavel[state.l_bins[k]],state.wavel[state.r_bins[k]]]),$
                0.0,ssi(k)
            endfor
            
            
            plot,state.wavel,fopt*skymask,xrange=xr,$
              title='optimal',yrange=[bmep_percent_cut(fopt,2),bmep_percent_cut(fopt,98)],/nodata
            oplot,state.wavel,fopt,color=254,psym=3
            index=where(fopt*skymask ne 0)
            oplot,state.wavel,fopt*skymask
            if overploterr then oploterr,state.wavel,fopt,fopterr,3
            for k=0,NBINS-1 do begin
              oplot,[state.wavel[state.l_bins[k]],state.wavel[state.l_bins[k]]],minmax(fopt),color=255
              oplot,[state.wavel[state.r_bins[k]],state.wavel[state.r_bins[k]]],minmax(fopt),color=255
              xyouts,avg([state.wavel[state.l_bins[k]],state.wavel[state.r_bins[k]]]),$
                0.0,ssi(k)
            endfor
            
            ;show the sky mask
            ;                if showskymask then begin
            ;                  maskedimg=state.data
            ;                  for jj=0,n_elements(maskedimg[*,0])-1 do maskedimg[jj,*]=maskedimg[jj,*]*skymask[jj]
            ;                  newimg=bytscl(maskedimg,top=255,/NAN,$
            ;                     min=bmep_percent_cut(maskedimg,10),max=bmep_percent_cut(maskedimg,90))
            ;                  if n_elements(im) gt 0 then im.setdata.newimg $
            ;                    else im=image(newimg,dimensions=[n_elements(newimg[*,0]),n_elements(newimg[0,*])],$
            ;                        window_title=slitname,margin=[0,0,0,0],max_value=255,min_value=0)
            ;                  endif
            
            ;find statistics of a certain area
            if leftstatspos ne -1.0 then begin
              bmep_plot_stats,state,leftstatspos,rightstatspos,f,fopt,statshift,ferr,fopterr
              index=where(state.wavel ge leftstatspos and state.wavel le rightstatspos,ct)
              if ct le 3 then print,'no flux between cursor marks' else begin
                if findautowidth eq 1 then begin
                  bmep_auto_width_calculator,j,centerarr,state,order,bkgndl,bkgndr,$
                    printp,fitgaussp,plotp,slidep,pwindowsize,singlep,$
                    cosmic_sigma,n_iterate_cosmic,bkgnd_naverage,max_rays_per_col,$
                    F,ferr,Fopt,fopterr,img2d_nobkgnd,parr,sky_reslut_arr,sky_residuals,$
                    widtharr,leftstatspos,rightstatspos
                endif ; find auto width
              endelse ; ct gt 3
            endif;stats
            
            ;if it should save the spectrum and extraction parameters
            if saveflag then begin
            
              flux=f
              err=ferr
              flux_opt=fopt
              erropt=fopterr
              wavel=state.wavel
              slitname=state.slitname
              img2d=state.data
              flux_opt_masked=fopt*skymask
              flux_masked=f*skymask
              
              if n_elements(centerarr) eq 1 then suffix='' else suffix='.'+ssi(j+1)
              if objnum ne -1 then suffix='.'+ssi(objnum)

              ;                  if state.savetext eq 0 then begin

;              forprint,wavel,flux,err,$
;                textout=state.savepath+slitname+suffix+'.txt',comment="# wavelength(A) Flux Flux_error"
;              forprint,wavel,flux_opt,erropt,$
;                textout=state.savepath+slitname+suffix+'_optimal.txt',comment="# wavelength(A) Flux Flux_error"

                
              if n_elements(centerarr) eq 1 then suffix='' else suffix='.'+ssi(j+1)
              if j eq 0 then suffix=''
              if objnum ne -1 then suffix='.'+ssi(objnum)
              if objnum ne -1 and n_elements(centerarr) ne 1 then $
                print,'WARNING, EXTRACTING MULTIPLE THINGS WITH OBJNUM NOT -1 IS BAD'
              extrainfo1=state.extrainfo1
              extrainfo2=state.extrainfo2
              extrainfo3=state.extrainfo3
              

              ;extract original name from fits header.
              index=WHERE(extrainfo1 eq 'MSKNM',ct)
              if ct eq 0 then message,'ERROR: MASKNAME (MSKNM) IS NOT DEFINED IN EXTRAINFO1'
              fitsfilenm=extrainfo2[index]+'.'
              index=WHERE(extrainfo1 eq 'FILTNM',ct)
              if ct eq 0 then message,'ERROR: FILTERNAME (FILTNM) IS NOT DEFINED IN EXTRAINFO1'
              fitsfilenm=fitsfilenm+extrainfo2[index]+'.'
              index=WHERE(extrainfo1 eq 'SLITNM',ct)
              if ct eq 0 then message,'ERROR: SLITNAME (SLITNM) IS NOT DEFINED IN EXTRAINFO1'
              fitsfilenm=fitsfilenm+extrainfo2[index]+suffix
              fitsfilenm=fitsfilenm[0]
              
              
              FORWARD_FUNCTION bmep_make_hdr
              
              ;make hdr
              ;exten 0
              header=bmep_make_hdr('',extrainfo1,extrainfo2,extrainfo3,j,centerarr,widtharr,$
                gausschiarr,fitgaussp,cosmic_sigma,order,state,objnum,wavel,autoextractflag, $
                wbyhand,cbyhand,mom1,mom2,slitloss,$
                gwidth, gcenter, gamp, glinear, gwidth_err, gcenter_err, $
                gamp_err, glinear_err,bkgndl,bkgndr, /exten)
              writefits,state.savepath+fitsfilenm+'.1d.fits','',header
              
              ;exten 1
              header=bmep_make_hdr(flux_opt,extrainfo1,extrainfo2,extrainfo3,j,centerarr,widtharr,$
                gausschiarr,fitgaussp,cosmic_sigma,order,state,objnum,wavel,autoextractflag, $
                wbyhand,cbyhand,mom1,mom2,slitloss,$
                gwidth, gcenter, gamp, glinear, gwidth_err, gcenter_err, $
                gamp_err, glinear_err,bkgndl,bkgndr,/image)
              writefits,state.savepath+fitsfilenm+'.1d.fits',flux_opt,header,/append
              
              ;exten 2
              header=bmep_make_hdr(erropt,extrainfo1,extrainfo2,extrainfo3,j,centerarr,widtharr,$
                gausschiarr,fitgaussp,cosmic_sigma,order,state,objnum,wavel,autoextractflag, $
                wbyhand,cbyhand,mom1,mom2,slitloss,$
                gwidth, gcenter, gamp, glinear, gwidth_err, gcenter_err, $
                gamp_err, glinear_err,bkgndl,bkgndr,/image)
              writefits,state.savepath+fitsfilenm+'.1d.fits',erropt,header,/append
              
              ;exten 3
              header=bmep_make_hdr(flux,extrainfo1,extrainfo2,extrainfo3,j,centerarr,widtharr,$
                gausschiarr,fitgaussp,cosmic_sigma,order,state,objnum,wavel,autoextractflag, $
                wbyhand,cbyhand,mom1,mom2,slitloss,$
                gwidth, gcenter, gamp, glinear, gwidth_err, gcenter_err, $
                gamp_err, glinear_err,bkgndl,bkgndr,/image)
              writefits,state.savepath+fitsfilenm+'.1d.fits',flux,header,/append
              
              ;exten 4
              header=bmep_make_hdr(err,extrainfo1,extrainfo2,extrainfo3,j,centerarr,widtharr,$
                gausschiarr,fitgaussp,cosmic_sigma,order,state,objnum,wavel,autoextractflag, $
                wbyhand,cbyhand,mom1,mom2,slitloss,$
                gwidth, gcenter, gamp, glinear, gwidth_err, gcenter_err, $
                gamp_err, glinear_err,bkgndl,bkgndr,/image)
              writefits,state.savepath+fitsfilenm+'.1d.fits',err,header,/append
              
              ;exten 5
              header=bmep_make_hdr(state.ydata,extrainfo1,extrainfo2,extrainfo3,j,centerarr,widtharr,$
                gausschiarr,fitgaussp,cosmic_sigma,order,state,objnum,wavel,autoextractflag, $
                wbyhand,cbyhand,mom1,mom2,slitloss,$
                gwidth, gcenter, gamp, glinear, gwidth_err, gcenter_err, $
                gamp_err, glinear_err,bkgndl,bkgndr,/image,/no_wave_info)
              writefits,state.savepath+fitsfilenm+'.1d.fits',state.ydata,header,/append
              
              ;exten 6
              header=bmep_make_hdr(state.ydataerr,extrainfo1,extrainfo2,extrainfo3,j,centerarr,widtharr,$
                gausschiarr,fitgaussp,cosmic_sigma,order,state,objnum,wavel,autoextractflag, $
                wbyhand,cbyhand,mom1,mom2,slitloss,$
                gwidth, gcenter, gamp, glinear, gwidth_err, gcenter_err, $
                gamp_err, glinear_err,bkgndl,bkgndr,/image,/no_wave_info)
                
              writefits,state.savepath+fitsfilenm+'.1d.fits',state.ydataerr,header,/append
              revisevar=0
              print,'saved ',fitsfilenm
              
              
              
;                      outpath=state.savepath+fitsfilenm+'.1d.txt'
;                      forprint,wavel,flux,err,$
;                      textout=outpath,comment="# wavelength(A) Flux Flux_error"
;
;                      outpath=state.savepath+fitsfilenm+'.1d_optimal.txt'
;                      forprint,wavel,flux_opt,erropt,$
;                      textout=outpath,comment="# wavelength(A) Flux Flux_error"
              
              ;if object is a star, update the starlist.txt
              index=WHERE(extrainfo1 eq 'ISSTAR',ct)
              IF ct eq 1 and extrainfo2[index[0]] eq 1 and objnum eq -1 and j eq 0 then begin
                readcol,state.savepath+'00_starinfo.txt',maskstar,$
                  filtstar,objstar,yexpect_star,yactual_star,widthstar,sigmastar,/silent,format='A,A,A,F,F,F,F'
                
                index=WHERE(extrainfo1 eq 'MSKNM')
                themask=extrainfo2[index[0]]
                index=WHERE(extrainfo1 eq 'FILTNM')
                thefilter=extrainfo2[index[0]]
                index=WHERE(extrainfo1 eq 'SLITNM')
                theobj=extrainfo2[index[0]]
                index=where(extrainfo1 eq 'YEXPECT')
                theyexpect=extrainfo2[index[0]]
                thesigmastar=gwidth[j]
                
                index=where((maskstar eq themask) and (filtstar eq thefilter) $
                  and (objstar eq theobj),ct)
                if ct eq 1 then begin
                
                  yexpect_star[index]=theyexpect
                  yactual_star[index]=gcenter[j]
                  widthstar[index]=thesigmastar*2.0*SQRT(2.0*ALOG(2.0))
                  sigmastar[index]=thesigmastar
                endif else begin
                  maskstar=[maskstar,themask]
                  filtstar=[filtstar,thefilter]
                  objstar=[objstar,theobj]
                  yexpect_star=[yexpect_star,theyexpect]
                  yactual_star=[yactual_star,gcenter[j]]
                  widthstar=[widthstar,thesigmastar*2.0*SQRT(2.0*ALOG(2.0))]
                  sigmastar=[sigmastar,thesigmastar]
                endelse ; ct
                forprint,textout=state.savepath+'00_starinfo.txt',maskstar+' ',$
                  filtstar+' ',objstar+' ',yexpect_star,yactual_star,widthstar,sigmastar,/nocomment
                  
                  
              endif ; extrainfo2[index]
            ;                endif ; state.savetext (REMOVED)
            endif ; save flag
          endfor ; centerarr
          
          ;automatically turn off auto width.
          if  findautowidth eq 1 and leftstatspos ne -1.0 then findautowidth=0
          
        endif ; extraction
        ;automatically turn off saveing
        saveflag=0
        
        ;get user input
        key=get_kbrd()
        ;========================================================================================
        ;========================================================================================
        ;========================================================================================
        ;========================================================================================
        ;========================================================================================
        ;========================================================================================
        ;
        
        
        
        ;below this point are all of the keyboard commands...
        
        ;stop command for debugging purposes.
        if key eq '0' then begin
          stop
        endif
        
        ;1-5 change object number.
        if key eq '1' then begin
          objnum=-1
          print,'object number, ', objnum
        endif
        
        if key eq '2' then begin
          objnum=2
          print,'object number, ', objnum
        endif
        
        if key eq '3' then begin
          objnum=3
          print,'object number, ', objnum
        endif
        
        if key eq '4' then begin
          objnum=4
          print,'object number, ', objnum
        endif
        
        if key eq '5' then begin
          objnum=5
          print,'object number, ', objnum
        endif
        ;check objnum and multiple object conflict.
        if objnum ne -1 and n_elements(centerarr) gt 1 then $
          print,'please only one at a time if changing the object number'
          
        ;edit bkgnd_avg
        if key eq 'a' then begin
          print,'current bkgnd_naverage ',bkgnd_naverage
          print,'enter new bkgnd_naverage'
          read,bkgnd_naverage
        endif
        
        ;toggle find auto width. SNR as function of width...
        if key eq 'A' then begin
          if findautowidth eq 1 then  findautowidth=0 $
          else if findautowidth eq 0 then findautowidth=1 $
        else findautowidth=0
      endif
      
      ;edit cosmic sigma threshold.
      if key eq 'c' then begin
        print,'current cosmic sigma ',cosmic_sigma
        print,'enter new cosmic sigma'
        read,cosmic_sigma
      endif
      
      ;edit max num rays per col
      if key eq 'C' then begin
        print,'current max_rays_per_col ',max_rays_per_col
        print,'enter new max_rays_per_col'
        read,max_rays_per_col
      endif
      
      ;delete center array and those associated w/ it
      IF key eq 'd' then begin
        bkgndl=[]
        bkgndr=[]
        centerarr=[]
        gausschiarr=[]
        wbyhand=[]
        cbyhand=[]
        widtharr=[]
        mom1=[]
        mom2=[]
        slitloss=[]
        gwidth=[]
        gcenter=[]
        gamp=[]
        glinear=[]
        gwidth_err=[]
        gcenter_err=[]
        gamp_err=[]
        glinear_err=[]
        print,'center and other arrays deleted'
      endif
      
      ;toggle extraction.
      if key eq 'e' then begin
        if extraction eq 1 then  extraction=0 $
        else if extraction eq 0 then extraction=1 $
      else extraction=0
    endif
    
    ;change sky mask percent
    if key eq 'f' then begin
      print,skymaskpercent
      print,'enter new sky percent'
      read,skymaskpercent
    endif
    
    ;toggle showing the skymask
    if key eq 'F' then begin
      if showskymask eq 1 then showskymask=0 $
      else if extraction eq 0 then showskymask=1 $
    else showskymask=0
  endif
  
  ;toggle gaussian fit to get P
  if key eq 'g' then begin
    if fitgaussp eq 1 then fitgaussp=0 $
    else if fitgaussp eq 0 then fitgaussp=1 $
  else fitgaussp=0
  if fitgaussp eq 1 then print,'FITTING TO A GAUSSIAN.'
  if fitgaussp eq 0 then print,'not fitting to gaussian.'
endif

;print help sorted by functionality
if key eq '?' then begin
  print,'Help sorted by functionality'
  print,'---help options---'
  print,'? - Help display sorted by functionality'
  print,'h - Help display sorted alphabetically'
  print,'0 - stop command for debugging'
  print
  print,'---Plotting options---'
  print,'F - Toggle on viewing the sky mask'
  print,'p - Toggle plotting p (shown in red)'
  print,'r - Get statistics from an area (use only after extracting)'
  print,'v - Toggle plotting errors to data plot'
  print,'x - Change x range plotted'
  print
  print,'---Extraction options---'
  print,'e - Toggle extraction of the spectra'
  print,'o - Edit the ORDER of the background fit'
  print,'c - Change sigma thresh (def = 4.0)'
  print,'f - Edit sky mask percentage (both unmasked and masked are saved)'
  print,'g - Toggle p extraction is fit to gaussian'
  print,'P - Edit p parameters'
  print,'i - Toggle flat to not automaticall extract later'
  print,'Z - Guess redshift'
  print,'X - Fit redshift.'
  print
  print,'---Adding or subtracting object options---'
  print,'m - Mark an object automatically (and do sky region)'
  print,'M - Change width for "m" fit'
  print,'n - Mark an object where the cursor is'
  print,'d - Delete all extraction points'
  print,'w - Decrease width of extraction area'
  print,'W - Increase width of extraction area'
  print,'A - Toggle Auto-width calculator (needs statistics (r))'
  print,'1 - reset object num'
  print,'2 - set object num to 2'
  print,'3 - set object num to 3'
  print,'4 - set object num to 4'
  print,'5 - set object num to 5'
  print
  print,'---Sky options---'
  print,'s - Set the background (once on the left and once on the right)'
  print,'z - Delete sky area (all)'
  print,'S - Toggle save the data (also extract if not selected)'
  print,'a - Change the b_naverage parameter'
endif

;print help sorted alphabetically
if key eq 'h' then begin
  print,'Help sorted alphabetically'
  print,'? - Help display sorted by functionality'
  print,'0 - stop command for debugging'
  print,'1 - reset object num'
  print,'2 - set object num to 2'
  print,'3 - set object num to 3'
  print,'4 - set object num to 4'
  print,'5 - set object num to 5'
  print,'a - Change the b_naverage parameter'
  print,'A - Toggle Auto-width calculator (needs statistics (r))'
  print,'c - Change sigma thresh (def = 4.0)'
  print,'d - Delete all extraction points'
  print,'e - Toggle extraction of the spectra'
  print,'f - Edit sky mask percentage'
  print,'F - Toggle on viewing the sky mask'
  print,'g - Toggle p extraction is fit to gaussian'
  print,'h - Help display sorted alphabetically'
  print,'i - Toggle flat to not automaticall extract later'
  print,'m - Mark an object automatically (and do sky region)'
  print,'M - Change width for "m" fit'
  print,'n - Mark an object where the cursor is'
  print,'o - Edit the ORDER of the background fit'
  print,'P - Edit with p parameters'
  print,'p - Toggle plotting p (shown in red)'
  print,'r - Get statistics from an area (use only after extracting)'
  print,'s - Set the background (once on the left and once on the right)'
  print,'S - Toggle save the data (also extract if not selected)'
  print,'v - Toggle add errors to data plot'
  print,'w - Decrease width of extraction area'
  print,'W - Increase width of extraction area'
  print,'x - Change x range plotted'
  print,'X - Fit redshift.'
  print,'z - Delete sky area (all)'
  print,'Z - Guess redshift'
endif

;change the automatically extract flag.
if key eq 'i' then begin
  if autoextractflag eq 1 then BEGIN
    autoextractflag=0
    print,'set to NOT automatically extract.'
  endif else begin
    autoextractflag=1
    print,'set to automatically extract.'
  endelse
endif

;change the mwidth
if key eq 'M' then begin
  print,mwidth
  print,'enter new mwidth:'
  read,mwidth
  if mwidth lt 0 then print,'no negative orders plz'
  mwidth=abs(mwidth)
  if round(mwidth) ne mwidth then print,'why isnt the mwidth an integer'
  mwidth=round(mwidth)
endif

;fit a gaussian and center on profile
if key eq 'm' then begin
  cursor,xpos,ypos,/nowait
  if xpos lt 0 then xpos = n_elements(state.ydata)/2
  
  lower=round(xpos-mwidth)
  if lower lt 0 then lower = 0
  upper=round(xpos+mwidth)
  if upper ge n_elements(state.ydata) then upper = n_elements(state.ydata)-1
  yfit=state.ydata[lower:upper]
  yfiterr=state.ydataerr[lower:upper]; error is SIGMA, not VARAIANCE.
  nterms=4
  
  
  pi =[{fixed:0, limited:[1,1], limits:[0.D,max(yfit)*1.2]},$ ;peak value
    {fixed:0, limited:[1,1], limits:[0.D,n_elements(state.ydata)]},$ ;peak centroid
    {fixed:0, limited:[0,0], limits:[0.D,0.D]},$ ;sigma
    {fixed:1, limited:[0,0], limits:[0.D,0.D]}];,$ ;linear bkgnd term
  ;{fixed:0, limited:[0,0], limits:[0.D,0.D]}]  ;quadratic background term
  estimates=[yfit[fix(xpos)-lower],xpos,2.0,0.0]
  oplot,findgen(upper-lower+1)+lower,yfit,color=255.+255.*255.
  xfit=findgen(upper-lower+1)+lower
  dummy=MPFITPEAK(xfit,yfit,$
    coeff,nterms=nterms,error=yfiterr,sigma=gauss_sigma,/gaussian,$
    estimates=estimates,parinfo=pi,status=status,chisq=chisq)
  print,'fit status (1 is good) ' ,status
  if status eq 1 or status eq 3 then begin
    oplot,findgen(upper-lower)+lower,dummy,color=245
    wait,1
    
    print,'peak    : ',coeff[0],' pm ',gauss_sigma[0],' estimate was ',estimates[0]
    print,'center  : ',coeff[1],' pm ',gauss_sigma[1],' estimate was ',estimates[1]
    print,'width   : ',coeff[2],' pm ',gauss_sigma[2],' estimate was ',estimates[2]
    print,'linear  : ',coeff[3],' pm ',gauss_sigma[3],' estimate was ',estimates[3]
    ;print,'quad    : ',coeff[4],' pm ',gauss_sigma[4],' estimate was ',estimates[4]
    fwhm_fit=2.0*SQRT(2.0*ALOG(2.0))*coeff[2]
    print,'2xfwhm   : ',fwhm_fit
    print,'r-chisqr : ',chisq/float(n_elements(yfit)-1.0)
    if coeff[1] gt 0 and coeff[1] lt n_elements(state.ydata) then begin
      ; fix rounding issues
      fwhm_fit=round(fwhm_fit)
      
      ;check min width
      index=where(state.extrainfo1 eq 'MINW',ct)
      if ct eq 1 then begin
        minw=state.extrainfo2[index[0]]
        print,'minw is ',minw
        print,'fwhmx2 is ',fwhm_fit
        if minw gt 0 and fwhm_fit lt minw then fwhm_fit=minw
        print,'fwhmx2 is ',fwhm_fit,' after minw correction'
      endif
      
      ;moments
      mom1=[mom1,total(xfit*yfit)/total(yfit)]
      mom2=[mom2,2.355*sqrt(abs(total(xfit*xfit*yfit)/total(yfit) - $
        (total(xfit*yfit)/total(yfit))*(total(xfit*yfit)/total(yfit))))]
      print,'moment 1 is ',mom1[-1]
      print,'moment 2 is ',mom2[-1]
      chisq=chisq/float(n_elements(yfit)-1.0)
      
      
      centerarr=fix([centerarr,round(coeff[1])])
      widtharr=round([widtharr,fwhm_fit])
      gausschiarr=[gausschiarr,chisq]
      wbyhand=[wbyhand,0]
      cbyhand=[cbyhand,0]
      slitloss=[slitloss,bmep_calc_slitloss(coeff[2])]
      gwidth=[gwidth,coeff[2]]
      gcenter=[gcenter,coeff[1]]
      print,gcenter
      gamp=[gamp,coeff[0]]
      glinear=[glinear,coeff[3]]
      gwidth_err=[gwidth_err,gauss_sigma[2]]
      gcenter_err=[gcenter_err,gauss_sigma[1]]
      gamp_err=[gamp_err,gauss_sigma[0]]
      glinear_err=[glinear_err,gauss_sigma[3]]
      
      ;if automatically determining the sky...
      if autosky then begin
        bkgndl=[coeff[1]-3.0*SQRT(2.0*ALOG(2.0))*coeff[2]-autoskybuffer-autoskywidth]
        bkgndr=[coeff[1]-3.0*SQRT(2.0*ALOG(2.0))*coeff[2]-autoskybuffer]
        bkgndl=[bkgndl,coeff[1]+3.0*SQRT(2.0*ALOG(2.0))*coeff[2]+autoskybuffer]
        bkgndr=[bkgndr,coeff[1]+3.0*SQRT(2.0*ALOG(2.0))*coeff[2]+autoskybuffer+autoskywidth]
        
        bkgndl=fix(bkgndl)
        bkgndr=fix(bkgndr)
        
      endif;auto sky
      
    endif ;coeff makes sense
  endif ;status eq 1
endif ; key eq m

;create a profile FORCING THE CENTER TO BE ON THE MOUSE
if key eq 'n' then begin
  cursor,xpos,ypos,/nowait
  if xpos gt 0 and xpos lt n_elements(state.ydata) then begin
    centerarr=[centerarr,fix(xpos)]
    widtharr=[widtharr,fix(default_width)]
    gausschiarr=[gausschiarr,-1.0]
    wbyhand=[wbyhand,1]
    cbyhand=[cbyhand,1]
    mom1=[mom1,-1.0]
    mom2=[mom2,-1.0]
    slitloss=[slitloss,bmep_calc_slitloss(default_width/2.355)]
    gwidth=[gwidth,-99]
    gcenter=[gcenter,-99]
    gamp=[gamp,-99]
    glinear=[glinear,-99]
    gwidth_err=[gwidth_err,-99]
    gcenter_err=[gcenter_err,-99]
    gamp_err=[gamp_err,-99]
    glinear_err=[glinear_err,-99]
    print,'centers ',centerarr
    print,'widths ',widtharr
  endif
endif

;change the order
if key eq 'o' then begin
  print,order
  print,'enter new order for bkgnd fit:'
  read,order
  if order lt 0 then print,'no negative orders plz'
  order=abs(order)
  if fix(order) ne order then print,'why isnt the order an integer'
  order=fix(order)
  if order gt 10 then print, 'really? why would you need an order that high.'
endif

;Allow user to change P parameters.
if key eq 'P' then begin
  pspot1:
  print,'use one value for p across all wavelengths? (set to: '+ssi(singlep)+') (type 1 or 0)'
  read,singlep
  
  if singlep eq 0 then begin
    slidep=1
    print,'set window size to bin in wavelength direction (set to: '+ssf(pwindowsize)+')
    read,pwindowsize
  endif
  
  if singlep eq 1 then begin
    slidep=0
    pspot2:
    print,'fit p to a gaussian? (set to:'+ssi(fitgaussp)+') (type 1 or 0)'
    read,fitgaussp
    if fitgaussp ne 1 and fitgaussp ne 0 then goto, pspot2
  endif
  
  if singlep ne 1 and singlep ne 0 then goto, pspot1
  
endif ; end of p options

;Toggle plotp and printp
if key eq 'p' then begin
  if printp eq 1 then  printp=0 $
  else if printp eq 0 then printp=1 $
else printp=0
if plotp eq 1 then  plotp=0 $
else if plotp eq 0 then plotp=1 $
else plotp=0
endif

;Get left and right spots for statics calculations
if key eq 'r' then begin
  cursor,leftstatspos,ypos,/nowait
  print,leftstatspos
  print,'hit r again'
  key=' '
  while key ne 'r' do key=get_kbrd()
  cursor,rightstatspos,ypos,/nowait
  print,rightstatspos
  ;check if backwards
  if leftstatspos gt rightstatspos then begin
    temp=leftstatspos
    leftstatspos=rightstatspos
    rightstatspos=temp
  endif ; backwards
endif;key r

;set the background
if key eq 's' then begin
  print,'hit s again'
  cursor,xleft,ypos,/nowait
  key=get_kbrd()
  while key ne 's' do key=get_kbrd()
  cursor,xright,ypos,/nowait
  if fix(xleft) ne fix(xright) then begin
    if xleft gt xright then begin
      tmp=xleft
      xleft=xright
      xright=tmp
    endif ; need to swap
    bkgndl=[bkgndl,xleft]
    bkgndr=[bkgndr,xright]
  endif ; xl and xr not eq
  print,'OK...'
endif ; key s



;save the spectrum
if key eq 'S' then begin
  if saveflag then saveflag=0 else begin
    print,'Saving the data.'
    saveflag=1
  endelse
  if extraction eq 0 then extraction=1
endif

;toggle plotting the error bars.
if key eq 'v' then begin
  if overploterr eq 1 then  overploterr=0 $
  else if overploterr eq 0 then overploterr=1 $
else overploterr=0
endif

;increase width of the extraction
if key eq 'w' then begin
  if n_elements(widtharr) gt 1 then begin
    print,'current widths: '
    print,widtharr
    print,'enter new widths:'
    for ii=0,n_elements(widtharr)-1 do begin
      temp=1
      read,temp,prompt=ssi(ii+1)+" >"
      if widtharr[ii] ne temp then wbyhand[ii]=1
      if widtharr[ii] ne temp then slitloss[ii]=bmep_calc_slitloss(temp/2.355)
      widtharr[ii]=temp
    endfor
  endif else begin
    widtharr=widtharr+1.0
    wbyhand[0]=1
    slitloss[0]=bmep_calc_slitloss(widtharr[0]/2.355)
  endelse
  print,'width: ', widtharr
endif

;Decrease width of the extraction
if key eq 'W' then begin
  if n_elements(widtharr) gt 1 then begin
    print,'current widths: '
    print,widtharr
    print,'enter new widths:'
    for ii=0,n_elements(widtharr)-1 do begin
      temp=1
      read,temp,prompt=ssi(ii+1)+" >"
      if widtharr[ii] ne temp then wbyhand[ii]=1 ; check if kept same value
      if widtharr[ii] ne temp then slitloss[ii]=bmep_calc_slitloss(temp/2.355)
      widtharr[ii]=temp
    endfor
  endif else begin
    widtharr=widtharr-1.0
    wbyhand[0]=1
    slitloss[0]=bmep_calc_slitloss(widtharr[0]/2.355)
  endelse
  ;check for extraction widths that are too small
  index=where(widtharr lt 1,ct)
  if ct gt 0 then begin
    widtharr[index]=1.0
    slitloss[index]=bmep_calc_slitloss(widtharr[index]/2.355)
    endif
  print,'width: ', widtharr
endif

;adjust x range of plot
if key eq 'x' then begin
  cursor,leftxpos,ypos,/nowait
  print,leftxpos
  print,'hit x again'
  key=' '
  while key ne 'x' do key=get_kbrd()
  cursor,rightxpos,ypos,/nowait
  print,rightxpos
  ;check if want to restart range
  if leftxpos ge rightxpos then begin
    leftxpos=min(state.wavel)
    rightxpos=max(state.wavel)
  endif ; backwards
  xr=[leftxpos,rightxpos]
;          temp=1.123
;          print,'enter xmin ',xr[0]
;          read,temp
;          xr[0]=temp
;          print,'enter xmax ',xr[1]
;          read,temp
;          xr[1]=temp
endif

;guess redshift for header purposes
if key eq 'X' then begin
  nofit=0
  !P.MULTI=[0,1,1]
  xr_save=xr
  cursor,leftxpos,ypos,/nowait
  print,leftxpos
  print,'hit X again'
  key=' '
  while key ne 'X' do key=get_kbrd()
  cursor,rightxpos,ypos,/nowait
  print,rightxpos
  ;check if want to restart range
  if leftxpos ge rightxpos then begin
    if leftxpos eq rightxpos then nofit=1 $
    else goto, Xexit ; GOTO STATEMENT, DEAL WITH IT.
  endif ; backwards
  xr=[leftxpos,rightxpos]
  
  ;wl of best center guess
  xpos = state.wavel[nearest_to(state.wavel,avg(xr))]
  
  lower=round(nearest_to(state.wavel,leftxpos))
  upper=round(nearest_to(state.wavel,rightxpos))
  yfit=fopt[lower:upper]
  yfiterr=fopterr[lower:upper]; error is SIGMA, not VARAIANCE.
  yfiterr=yfiterr/max(yfit)
  yfit=yfit/max(yfit)
  xfit=state.wavel[lower:upper]
  nterms=4
  
  plot,xfit,yfit
  
  
  pi =[{fixed:0, limited:[1,1], limits:[min(yfit),max(yfit)*1.2]},$ ;peak value
    {fixed:0, limited:[1,1], limits:[double(xr)]},$ ;peak centroid
    {fixed:0, limited:[0,0], limits:[0.D,0.D]},$ ;sigma
    {fixed:0, limited:[0,0], limits:[0.D,0.D]}];,$ ;linear bkgnd term
  ;{fixed:0, limited:[0,0], limits:[0.D,0.D]}]  ;quadratic background term
  estimates=[max(yfit),xpos,2.0,min(yfit)]
  oplot,xfit,yfit,color=255.+255.*255.
  dummy=MPFITPEAK(xfit,yfit,$
    coeff,nterms=nterms,error=yfiterr,sigma=gauss_sigma,/gaussian,$
    estimates=estimates,parinfo=pi,status=status,chisq=chisq)
  print,'fit status (1 is good) ' ,status
  if status eq 1 or status eq 3 then begin
    oplot,xfit,dummy,color=245
    oplot,[coeff[1],coeff[1]],minmax(yfit)
    
    
    print,'peak    : ',coeff[0],' pm ',gauss_sigma[0],' estimate was ',estimates[0]
    print,'center  : ',coeff[1],' pm ',gauss_sigma[1],' estimate was ',estimates[1]
    print,'width   : ',coeff[2],' pm ',gauss_sigma[2],' estimate was ',estimates[2]
    print,'linear  : ',coeff[3],' pm ',gauss_sigma[3],' estimate was ',estimates[3]
    ;print,'quad    : ',coeff[4],' pm ',gauss_sigma[4],' estimate was ',estimates[4]
    
    
    ;define line names
    linenames= ['Ha',    'NII', '[OII]','[OIII]','[OIII]2','Hb',    'SII','hgamma']
;    linewavels=[6564.614,6585.27,3728.48,5008.239,4960.295,4862.721,6718.29,4339.00];vacuum
    linewavels=[6562.801,6583.45,3727.28,5006.843,4958.911,4861.363,6716.44,4341.00];air
    
    forprint,indgen(n_elements(linenames)),' '+linenames,linewavels
    print,n_elements(linenames),'other'
    choice=-1
    print,'enter choice'
    read,choice
    
    if choice ge 0 and choice le n_elements(linenames) then begin
      if choice eq n_elements(linenames) then begin
        new_name='?'
        new_wavel=-99.99
        print,'enter line name
        read,new_name
        print,'enter line wavelength (make notice of if u need vac/air)'
        read,new_wavel
      endif else begin
        new_name=linenames[choice]
        new_wavel=linewavels[choice]
      endelse
      
      redshift=(coeff[1]/new_wavel)-1.0
      ;redshifterr=redshift*(gauss_sigma[1]/new_wavel)
      ;                redshifterr=avg([,((coeff[1]-gauss_sigma[1])/new_wavel)-1.0])
      redshifterr=abs(((coeff[1]+gauss_sigma[1])/new_wavel)-1.0-redshift)

      print,'z=',redshift
      print,'line name=',new_name
      print,'added to fits hdr.'
      
      if objnum ne -1 then slitname=state.slitname+'---'+ssi(objnum) else slitname=state.slitname
      index=where(state.extrainfo1 eq 'MSKNM',ct)
      IF ct EQ 1 then maskname=state.extrainfo2[index[0]] else maskname='UNKNOWN'
      
      print,maskname
      print,'is the maskaname'
      
      if ~file_test('~/REDSHIFT_CATALOG_bmep.txt') then $
        forprint,maskname+' ',slitname,redshift,redshifterr,' '+new_name,new_wavel,coeff[1],$
        textout='~/REDSHIFT_CATALOG_bmep.txt',comment="# maskname slit z zerr linename restwave obswave", $
        format='(A20,A14,F10.6,F13.8,A12,F11.3,F11.3)' $
      else begin
      readcol,'~/REDSHIFT_CATALOG_bmep.txt',v1,v2,v3,v4,v5,v6,v7,format='A,A,F,F,A,F,F'
      index=where(sss(v1) eq sss(maskname) and sss(v2) eq sss(slitname) and sss(v5) eq sss(new_name), ct)
      ;                    forprint,v1,v2,v5
      ;                    print,maskname,slitname,new_name
      print,ct,' exact matches found..'
      if ct ge 1 then begin
        v1[index]=maskname
        v2[index]=slitname
        v3[index]=redshift
        v4[index]=redshifterr
        v5[index]=new_name
        v6[index]=new_wavel
        v7[index]=coeff[1]
      endif else begin
        v1=[v1,maskname]
        v2=[v2,slitname]
        v3=[v3,redshift]
        v4=[v4,redshifterr]
        v5=[v5,new_name]
        v6=[v6,new_wavel]
        v7=[v7,coeff[1]]
      endelse
      
      ;SORT THE LIST ->MASK -> SLIT
      index=bsort(v2)
      v1=v1[index]
      v2=v2[index]
      v3=v3[index]
      v4=v4[index]
      v5=v5[index]
      v6=v6[index]
      v7=v7[index]
      
      index=bsort(v1)
      v1=v1[index]
      v2=v2[index]
      v3=v3[index]
      v4=v4[index]
      v5=v5[index]
      v6=v6[index]
      v7=v7[index]
      
      forprint,v1+' ',v2,v3,v4,' '+v5,v6,v7,$
        textout='~/REDSHIFT_CATALOG_bmep.txt',comment="# maskname slit z zerr linename restwave obswave",$
        format='(A20,A14,F10.6,F13.8,A12,F11.3,F11.3)'
      index=where(sss(v1) eq sss(maskname) and sss(v2) eq sss(slitname),ct)
      if ct ge 1 then $
        forprint,v1[index]+' ',v2[index],v3[index],v4[index],' '+v5[index],v6[index],v7[index],$
          comment="# maskname slit z zerr linename restwave obswave",$
          format='(A20,A14,F10.6,F13.8,A12,F11.3,F11.3)'
    endelse ; file found
  endif;choice
endif ;status eq 1
xr=xr_save
Xexit:
endif;key X

;delete background arrays.
IF key EQ 'z' THEN BEGIN

  bkgndl=[]
  bkgndr=[]
  print,'background arrays deleted'
endif

;guess redshift
IF key EQ 'Z' THEN BEGIN

  ;            print,'there are ',state.n_bins,' number of bins'
  FORWARD_FUNCTION bmep_guess_redshift
  zguess=bmep_guess_redshift(state,wavel,fopt,linename=linename,linecol=linecol)
  
endif





endwhile ; key ne q or Q

print,'exited the background fitting, extraction and whatnot.'
end
'r': begin
  dummy=replicate(0.d,n_elements(state.ydata))
  state.ydata=dummy
  state.ydataerr=dummy
  state.n_bins=0
  state.l_bins=[0,0,0,0,0,0,0,0,0,0,0]
  state.r_bins=[0,0,0,0,0,0,0,0,0,0,0]
  print,'reset cross section'
end
'c': begin
  if state.cont_mode eq 0 then begin
    state.cont_mode=1
    print,'cont mode on'
  endif else begin
    state.cont_mode=0
    print,'cont mode off'
  endelse
end

's': begin
  if state.stats_mode eq 0 then begin
    state.stats_mode=1
    print,'cont mode on'
  endif else begin
    state.stats_mode=0
    print,'cont mode off'
  endelse
end

'h': begin
  print,'begin help'
  print,'b - go into extraction mode'
  print,'r - reset cross section'
  print,'c - continuum mode'
  print
end

ELSE: begin
  print,'begin help'
  print,'b - go into extraction mode'
  print,'r - reset cross section'
  print,'c - continuum mode'
  print
end
ENDCASE

!p.multi=[0,1,1]
oWin.uvalue=state
RETURN, 0 ; Skip default event handling
END ; keyboard handler

;laplacian of a gaussian function
function bmep_LoG,x,y,sigma
return,(-1.0/(!pi*sigma^4))*(1.0-(x^2+y^2)/(2.0*sigma^2))*exp((-(x^2+y^2))/(2.0*sigma^2))
end






pro bmep_lris,path_to_dropbox=path_to_dropbox,path_to_output=path_to_output
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  ;if not !textout then astrolib
  astrolib
  if not keyword_set(path_to_dropbox) then path_to_dropbox='~/Dropbox/'
  if not keyword_set(path_to_output) then path_to_output='~/LRIS/'
  cd,path_to_output,current=original_dir
  
  message,'did you add extrainfo1 with maskname, filtername, and slitname?'
  
  ;clear away all windows!!
  close,/all
  ireset,/no_prompt
  
  
  ;find folders of masks..
  spawn,'pwd'
  spawn,'ls > tempfiles.txt'
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt
  temp=[]
  for i=0,n_elements(filenames)-1 do $
    if bmep_dir_exist(filenames[i]) then temp=[temp,filenames[i]]
  filenames=temp
  
  forprint,indgen(n_elements(filenames)),replicate(' ',n_elements(filenames)),filenames
  print,'which number'
  read,choice
  if choice lt 0 then goto,theend
  if choice ge n_elements(filenames) then goto,theend
  
  
  ;foldername is same as the mask name
  maskname=filenames[choice]
  print,'folder: ',maskname
  cd,maskname
  
  savepath=path_to_output+maskname+'/1d_extracted/'
  if not bmep_DIR_EXIST(savepath) then spawn,'mkdir 1d_extracted'
  
  ;check if fits files exist
  ;lris files have the format cutout_slitname.fits
  spawn,'ls cutout_*.fits > xx__spec_list.txt'
  readcol,'xx__spec_list.txt',filenames,format='A',/silent
  spawn,'rm xx__spec_list.txt
  if n_elements(filenames) eq 0 then message,'no fits files found... probably no cutouts'
  
  
  ;get the names of the slits
  slitnames=[]
  for i=0,n_elements(filenames)-1 do slitnames=[slitnames,bmep_get_slitname(filenames[i],/lris)]
  
  
  
  ;get which slit to do.
  forprint,indgen(n_elements(slitnames)),replicate('. ',n_elements(slitnames)),slitnames
  print,n_elements(filenames),' all files'
  print,'which slit?'
  print
  choice=0
  read,choice
  if choice lt 0 then goto,theend
  if choice gt n_elements(filenames) then goto,theend
  if choice eq n_elements(filenames) then choicearr=indgen(n_elements(filenames)) else choicearr=[choice]
  
  for i=0,n_elements(choicearr)-1 do begin
  
  
    slitname=slitnames[choicearr[i]]
    
    ;check if the files exist
    scifile='cutout_'+slitname+'.fits'
    if ~file_test(scifile) then print,'file '+scifile+' does not exist'
    
    if file_test(scifile) and slitname ne '' then begin
    
      ;read in files
      sciimg=readfits(scifile,shdr, /SILENT)
      ny=n_elements(sciimg[0,*])
      nx=n_elements(sciimg[*,0])
      
      ;get wavel calib
      ;set gain and n_images...
      case maskname of
        'a_1689_feb2012_blue': begin
          readcol,path_to_dropbox+'siana_group/bill/a_1689_2012_feb/pix_to_wave_feb_2012.txt',wavel
          n_images=6.0
          gain=4.0
          filtername='blue'
        end
        'may_2010_version_1':  begin
          readcol,path_to_dropbox+'siana_group/bill/a_1689_2010_may/pix_to_wave_MAY_v1.txt',wavel
          n_images=10.0
          gain=4.0
          filtername='blue'
        end
        'macs0717_feb2012':  begin
          readcol,path_to_dropbox+'siana_group/bill/macs0717_feb2012/pix_to_wave_macs0717.txt',wavel
          n_images=7.0
          gain=4.0
          filtername='blue'
        end
        'may_2010_standard':  begin
          readcol,'/Users/bill/LRIS/may_2010_standard/wolf_1346_wavel.txt',wavel
          n_images=1.0
          gain=4.0
          filtername='blue'
        end
        'feb_2012_standard':  begin
          readcol,'/Users/bill/LRIS/feb_2012_standard/feb_2012_standard_wavecal.txt',wavel
          n_images=1.0
          gain=4.0
          filtername='blue'
        end
        'may_2010_red_t':  begin
          readcol,'/Users/bill/LRIS/may_2010_red_t/wavel_calib.txt',wavel
          n_images=19.0
          gain=1.2
          filtername='red'
        end
        'may_2010_red_b':  begin
          readcol,'/Users/bill/LRIS/may_2010_red_b/wavel_calib.txt',wavel
          n_images=19.0
          gain=1.2
          filtername='red'
        end
        'may_2010_red_standard':  begin
          readcol,'/Users/bill/LRIS/may_2010_red_standard/wavel_calib.txt',wavel
          n_images=1.0
          gain=1.2
          filtername='red'
        end
        'a_1689_feb2012_red':  begin
          readcol,'/Users/bill/LRIS/a_1689_feb2012_red/wavel_calib.txt',wavel
          n_images=15.0
          gain=1.2
          filtername='red'
        end
        'a_1689_feb2012_red_std':  begin
          readcol,'/Users/bill/LRIS/a_1689_feb2012_red_std/wavel_calib.txt',wavel
          n_images=1.0
          gain=1.2
          filtername='red'
        end
        else: begin
          print,'WARINING, NO MATCHING WAVELENGTH CALIBRATION...'
          print,'WARINING, NO MATCHING WAVELENGTH CALIBRATION...'
          print,'WARINING, NO MATCHING WAVELENGTH CALIBRATION...'
          
          wavel=findgen(nx)
          n_images=1.0
          gain=1.0
          filtername='none'
        END
      endcase
      if n_elements(wavel) lt nx then message,'bad wavelength calibration read in'
      
      ;convert image to photons.
      sciimg=sciimg*gain*n_images ;4.0 is the gain.
      
      ;calculate variance image
      var_img=sciimg*gain*n_images
      
      
      botpercent=10.0
      toppercent=90.0
      index=where(sciimg ne 0)
      highval=bmep_percent_cut(sciimg[index],toppercent)
      lowval=bmep_percent_cut(sciimg[index],botpercent)
      
      
      ;RUN THE BMEP!
      bmep_display_image,sciimg,sciimg,var_img,highval,lowval,slitname,filtername,wavel,savepath,revisevar=1
    endif ; if files exist
  endfor ; choicearr
  
  
  theend:
  
  cd,original_dir
  print,'end of best mosfire extraction program'
end




function bmep_make_hdr,data_make,extrainfo1,extrainfo2,extrainfo3,j,centerarr,widtharr,$
    gausschiarr,fitgaussp,cosmic_sigma,order,state,objnum,wavel,autoextractflag, $
    wbyhand,cbyhand,mom1,mom2,slitloss,$
    gwidth, gcenter, gamp, glinear, gwidth_err, gcenter_err, gamp_err, glinear_err,$
    bkgndl,bkgndr,exten=exten,image=image,no_wave_info=no_wave_info



  ;make fits header for new 1d file!
  MKHDR, header, data_make, exten=exten, image=image
  FOR jj=0,n_elements(extrainfo1) -1 do $
    if strnumber(extrainfo2[jj]) then $
    if float(extrainfo2[jj]) eq fix(extrainfo2[jj]) then $
    sxaddpar, Header,extrainfo1[jj],fix(extrainfo2[jj]),extrainfo3[jj] else $
    sxaddpar, Header,extrainfo1[jj],float(extrainfo2[jj]),extrainfo3[jj] else $
    sxaddpar, Header,extrainfo1[jj],STRCOMPRESS(extrainfo2[jj], /REMOVE_ALL),extrainfo3[jj]
  ;add in all parameters for the fit.
  sxaddpar, Header, 'YPOS', centerarr[j],' Y position of center of extraction'
  sxaddpar, Header, 'WIDTH', widtharr[j],' Half width in n_pixels'
  sxaddpar, Header, 'GAUSRCHI', gausschiarr[j],' Reduced Chi sqr of gauss fit to Y profile'
  sxaddpar, Header, 'WBYHAND', wbyhand[j],' Flag if width adjusted by hand'
  sxaddpar, Header, 'CBYHAND', cbyhand[j],' Flag if center adjusted by hand'
  sxaddpar, Header, 'MOM1', mom1[j],' First moment calculation'
  sxaddpar, Header, 'MOM2', mom2[j],' Second moment calc (transformed to FWHM)'
  sxaddpar, Header, 'SLITLOSS', slitloss[j], ' Flux in slit / Total flux assuming 2D Gaussians'
  sxaddpar, Header, 'GWIDTH', gwidth[j], ' Sigma of gaussian fit'
  sxaddpar, Header, 'GCENT', gcenter[j], ' Central pixel of gaussian fit'
  sxaddpar, Header, 'GAMP', gamp[j], ' Amplitude of gaussian fit'
  sxaddpar, Header, 'GLINEAR', glinear[j], ' Linear continuum of gaussian fit'
  sxaddpar, Header, 'GWIDTH_E', gwidth_err[j], ' 1sigma Error in Sigma of gaussian fit'
  sxaddpar, Header, 'GCENT_E', gcenter_err[j], ' 1sigma Error in CENTRAL pixel of gaussian fit'
  sxaddpar, Header, 'GAMP_E', gamp_err[j],  ' 1sigma Error in Amplitude of gaussian fit'
  sxaddpar, Header, 'GLINE_E', glinear_err[j], ' 1sigma Error in Linear continuum of gaussian fit'
  sxaddpar, Header, 'GAUSSP', fitgaussp,' Flag for if profile is fit to gaussian'
  sxaddpar, Header, 'CSIGMA', cosmic_sigma,' Sigma for cosmic ray rejection'
  sxaddpar, Header, 'ORDER', order,' Order of background fit (not for mosfire data)'
  sxaddpar, Header, 'NBKGND', order,' Number of areas that determine the background'
  FOR jj=0,n_elements(bkgndl)-1 do begin
    sxaddpar, Header, 'LBKGND'+ssi(jj+1), bkgndl[jj],' Left pixel number of background area'
    sxaddpar, Header, 'RBKGND'+ssi(jj+1), bkgndr[jj],' Right pixel number of background area'
  endfor
  sxaddpar, Header, 'NBINS', state.n_bins,' Number of areas of collapsed columns'
  FOR jj=0,state.n_bins-1 do begin
    sxaddpar, Header, 'LBIN'+ssi(jj+1), state.l_bins[jj],' Left pixel number of area of collapsed columns'
    sxaddpar, Header, 'RBIN'+ssi(jj+1), state.r_bins[jj],' Right pixel number of area of collapsed columns'
  endfor
  FOR jj=0,state.n_bins-1 do begin
    sxaddpar, Header, 'LBINW'+ssi(jj+1), wavel[state.l_bins[jj]],' Left wavelength of area of collapsed columns'
    sxaddpar, Header, 'RBINW'+ssi(jj+1), wavel[state.r_bins[jj]],' Right wavelength of area of collapsed columns'
  endfor
  
  sxaddpar, Header, 'CONTMODE', state.cont_mode,' Flag for "continuum mode"'
  if objnum eq -1 then sxaddpar, Header, 'OBJNUM', J+1, $
    ' Object number. 1=primary, >2 non-primary' $
  else sxaddpar, Header, 'OBJNUM', objnum,$
  ' Object number. 1=primary, >2 non-primary'
sxaddpar, Header, 'AUTOEX', autoextractflag,'Flag if should re-extract this object'
sxaddpar, Header, 'BLIND', 0, ' Flag if this extraction was done blindly'





sxaddpar, Header, 'COMMENT',' Exten 1: Optimal extraction'
sxaddpar, Header, 'COMMENT',' Exten 2: Optimal extraction error bars'
sxaddpar, Header, 'COMMENT',' Exten 3: Boxcar extraction'
sxaddpar, Header, 'COMMENT',' Exten 4: Boxcar extraction error bars'
sxaddpar, Header, 'COMMENT',' Exten 5: Light profile'
sxaddpar, Header, 'COMMENT',' Exten 6: Light profile error bars' 
sxaddpar, Header, 'COMMENT',' Anything called a "flag" means 1=yes and 0=no'

if keyword_set(no_wave_info) then begin
  sxaddpar,Header,'CRVAL1',1
  sxaddpar,Header,'CDELT1',1.0
  sxaddpar,Header,'CRPIX1',1
  sxaddpar,Header,'CTYPE1','LINEAR'
  endif


return,header


end



;make plots of the 1d extracted spectra
pro bmep_make_lris_plots,path_to_output=path_to_output
  astrolib
  
  ;defaults
  if not keyword_set(path_to_output) then path_to_output='~/LRIS/'
  cd,path_to_output,current=original_dir
  
  ;find folders of masks..
  spawn,'ls > tempfiles.txt'
  readcol,'tempfiles.txt',foldernames,format='A',/silent
  spawn,'rm tempfiles.txt'
  
  ;use iii as to not mess with restored variables
  for iii=0,n_elements(foldernames)-1 do begin
    extracted_folder=foldernames[iii]+'/1d_extracted'
    if bmep_dir_exist(extracted_folder) then begin
      cd,extracted_folder
      
      ;set filenames to nothing because if readcol
      ;does not find any files, the ones from the
      ;previous directory will still be in this var
      ;and it will crash.
      filenames=[]
      
      ;find .sav files
      spawn,'ls *.sav > tempfiles.txt'
      readcol,'tempfiles.txt',filenames,format='A',/silent,count=count
      spawn,'rm tempfiles.txt
      
      ;set flag to indicate if the PS file was started
      psstartflag=0
      
      ;loop through each save file
      for k=0,n_elements(filenames)-1 do begin
      
        ;begin PS file if this is first file for this folder.
        if k eq 0 then PS_Start, Filename='all_plots.ps'
        if k eq 0 then psstartflag=1
        
        print,'plotting... ',extracted_folder+filenames[k]
        if file_test(filenames[k]) then restore,filenames[k] else message,'no save file'
        
        ;if this is a good save file, it will have a variable called state.
        if n_elements(state) ne 0 then begin
          !p.multi=[0,1,2]
          cgplot,wavel,flux_opt,title=state.slitname,/xs,$
            ytitle='Optimal Flux',xtitle='Wavelength (Angstroms)'
          cgplot,wavel,flux,    title=state.slitname,/xs,$
            ytitle='Flux',xtitle='Wavelength (Angstroms)'
          !p.multi=[0,1,1]
        endif ; state ne 0
      endfor
      if psstartflag then begin
        ps_end
        psstartflag=0
        print,'ps ended'
      endif
    endif ; dir exists
    cd,path_to_output
  endfor ; filenames
end





;;read in redshift data
;readcol,path_to_output+'/allslits.seplines.zmosfire.qflag.txt',zmask,zslit,zobject,zredshift,format='X,A,I,I,F,X,X,X,X',/SILENT
;
;


pro bmep_mosdef_fluxcal

  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  !except=2 ;see division by zero errors instantly.
  astrolib
  
  
  
  ;creat list of 1d spec
  
  if not keyword_set(path_to_output) then path_to_output='~/mosfire/output/idl_output/2D' ; no trailing slash
  cd,path_to_output+'/1d_extracted',current=original_dir
  
  ;read in list of 1d spectra
  spawn,'rm *.fc.*
  spawn,'ls *.1d.fits > tempfiles.txt'
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt'
  print,'list contains ',n_elements(filenames),' objects
  
  ;read in normalization data
  readcol,path_to_output+'/1d_extracted/00_normalizations.txt',normmask,normfilt,normfactor,format='A,A,X,X,F',/SILENT
  
  ;loop through 1d spectra
  for i=0,n_elements(filenames)-1 do begin
    temp=readfits(filenames[i],hdr0,exten_no=0,/silent)
    ydata=readfits(filenames[i],hdr,exten_no=5,/silent)
    fopt=readfits(filenames[i],hdr,exten_no=1,/silent)
    fopterr=readfits(filenames[i],hdr,exten_no=2,/silent)
    f=readfits(filenames[i],hdr,exten_no=3,/silent)
    ferr=readfits(filenames[i],hdr,exten_no=4,/silent)
    
    ny=n_elements(ydata)
    
    
    maskname=strcompress(sxpar(hdr,'MSKNM'),/remove_all)
    filtername=strcompress(sxpar(hdr,'FILTNM'),/remove_all)
    slitname=strcompress(sxpar(hdr,'SLITNM'),/remove_all)
    
    
    INDEX=WHERE(normmask eq maskname and normfilt eq filtername,ct)
    
    if ct eq 1 then begin
      nfactor=normfactor[index[0]]
      fopt=fopt*nfactor
      fopterr=fopterr*nfactor
      f=f*nfactor
      ferr=ferr*nfactor
      
      sxaddpar,hdr0,'NFACTOR',nfactor,' Normalization factor'
      sxaddpar,hdr,'NFACTOR',nfactor,' Normalization factor'
      
      savename=strmid(filenames[i],0,strlen(filenames[i])-7)+'fc.1d.fits
      print,savename
      
      writefits,savename,'',hdr0
      writefits,savename,fopt,hdr,/append
      writefits,savename,fopterr,hdr,/append
      writefits,savename,f,hdr,/append
      writefits,savename,ferr,hdr,/append
      
      sxaddpar,hdr,'CRVAL1',1
      sxaddpar,hdr,'CDELT1',1.0
      sxaddpar,hdr,'CRPIX1',1
      sxaddpar,hdr,'CTYPE1','LINEAR'
      writefits,savename,ydata,hdr,/append
    endif;ct gt 0
  endfor ; filename
  cd,original_dir
end





pro bmep_mosdef_update_yexpect
  stop ;hammer time
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  !except=2 ;see division by zero errors instantly.
  astrolib
  
  ;creat list of 1d spec
  
  if not keyword_set(path_to_output) then path_to_output='~/mosfire/output/idl_output/2D' ; no trailing slash
  cd,path_to_output+'/1d_extracted',current=original_dir
  
  ;read in list of 1d spectra
  spawn,'ls *.1d.fits > tempfiles.txt'
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt'
  print,'list contains ',n_elements(filenames),' objects
  
  ;loop through 1d spectra
  for i=0,n_elements(filenames)-1 do begin
    temp=readfits(filenames[i],hdr0,exten_no=0,/silent)
    ydata=readfits(filenames[i],hdr,exten_no=5,/silent)
    fopt=readfits(filenames[i],hdr,exten_no=1,/silent)
    fopterr=readfits(filenames[i],hdr,exten_no=2,/silent)
    f=readfits(filenames[i],hdr,exten_no=3,/silent)
    ferr=readfits(filenames[i],hdr,exten_no=4,/silent)
    
    ny=n_elements(ydata)
    
    
    maskname=strcompress(sxpar(hdr,'MSKNM'),/remove_all)
    filtername=strcompress(sxpar(hdr,'FILTNM'),/remove_all)
    slitname=strcompress(sxpar(hdr,'SLITNM'),/remove_all)
    twodfilename=strcompress(sxpar(hdr,'FILNM'),/remove_all)
    
    slitlistfile=path_to_output+'/00mask_info/'+maskname+'_SlitList.txt'
    if ~file_test(slitlistfile) then message,'SLITLIST FILE MISSING... '+SLITLISTFILE
    
    ;1 2 17  45.61 -5  10  9.37  0.70  7.01  23763 800.00  -0.38 2 17  45.59 -5  10  9.65
    readcol,slitlistfile,slitnamearr,priorityarr,offsetarr,format='X,X,X,X,X,X,X,X,X,I,F,F,X,X,X,X,X,X',/silent
    index=where(slitnamearr eq slitname,ct)
    yexpect=-1
    isstar=-1
    minwidth=-1
    yshift=0.0
    if ct eq 1 then begin
      pixscale=0.1799
      midpoint=ny/2
      yexpect=midpoint+offsetarr[index[0]]/pixscale
      ;check if object is a star
      if abs(priorityarr[index]) eq 1 then isstar=1 else begin
        isstar=0
        readcol,path_to_output+'/1d_extracted/00_starinfo.txt',maskstar,$
          filtstar,objstar,yexpect_star,yactual_star,widthstar,sigmastar,/silent,format='A,A,A,F,F,F,F'
        index=where(maskstar eq maskname and filtstar eq filtername,ct)
        if ct ne 0 then begin
          print,ct,' number of stars found for ',maskname,' ',filtername
          yshift=avg(yexpect_star[index] - yactual_star[index])
          minwidth=fix(min(widthstar[index]))
        endif else print,'no matching stars found ',maskname,' ',filtername
      endelse ;if isstar
    endif
    
    yexpect=yexpect-yshift
    yexpect=round(yexpect)
    
    
    TEMP=readfits(path_to_output+'/'+twodfilename,twodhdr,exten_no=1,/silent)
    
    sxaddpar,hdr0,'CRVAL1',sxpar(twodhdr,'CRVAL1'),/savecomment
    sxaddpar,hdr0,'CDELT1',sxpar(twodhdr,'CDELT1'),/savecomment
    sxaddpar,hdr0,'CRPIX1',sxpar(twodhdr,'CRPIX1'),/savecomment
    sxaddpar,hdr0,'EXPTIME',sxpar(twodhdr,'EXPTIME'),/savecomment
    
    sxaddpar,hdr,'CRVAL1',sxpar(twodhdr,'CRVAL1'),/savecomment
    sxaddpar,hdr,'CDELT1',sxpar(twodhdr,'CDELT1'),/savecomment
    sxaddpar,hdr,'CRPIX1',sxpar(twodhdr,'CRPIX1'),/savecomment
    sxaddpar,hdr,'EXPTIME',sxpar(twodhdr,'EXPTIME'),/savecomment
    
    
    sxaddpar,hdr0,'YEXPECT',yexpect,/savecomment
    sxaddpar,hdr,'YEXPECT',yexpect,/savecomment
    sxaddpar,hdr0,'MINW',minwidth,/savecomment
    sxaddpar,hdr,'MINW',minwidth,/savecomment
    
    writefits,filenames[i],'',hdr0
    writefits,filenames[i],fopt,hdr,/append
    writefits,filenames[i],fopterr,hdr,/append
    writefits,filenames[i],f,hdr,/append
    writefits,filenames[i],ferr,hdr,/append
    
    sxaddpar,hdr,'CRVAL1',1
    sxaddpar,hdr,'CDELT1',1.0
    sxaddpar,hdr,'CRPIX1',1
    sxaddpar,hdr,'CTYPE1','LINEAR'
    writefits,filenames[i],ydata,hdr,/append
    
  endfor
  cd,original_dir
end





pro bmep_mosdef_rereduce
  stop
  ;setup
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  !except=2 ;see division by zero errors instantly.
  astrolib
  starttime=systime(/seconds)
  mwidth=8 ; width to fit a gaussian to...
  !p.multi=[0,1,1]
  ;spawn,'rm reex_*'
  
  ;creat list of 1d spec
  
  if not keyword_set(path_to_output) then path_to_output='~/mosfire/output/idl_output/2D' ; no trailing slash
  cd,path_to_output+'/1d_extracted',current=original_dir
  
  ;read in list of 1d spectra
  spawn,'ls *.1d.fits > tempfiles.txt' ;;for all files
;  spawn,'ls co3_01*.1d.fits > tempfiles.txt' ;for updating only one
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt'
  print,'list contains ',n_elements(filenames),' objects
  
  ;create backup folder
  foldermade=0
  prefixnum=0
  while foldermade eq 0 do begin
    if prefixnum lt 10 then $
      backup_folder_name='00'+ssi(prefixnum)+'backup' else $
      if prefixnum lt 100 then $
      backup_folder_name= '0'+ssi(prefixnum)+'backup' else $
      backup_folder_name=  ''+ssi(prefixnum)+'backup'
      
    if bmep_dir_exist(backup_folder_name) then print,backup_folder_name+' exists' else begin
      spawn,'mkdir '+backup_folder_name
      foldermade=1
    endelse
    prefixnum++
  endwhile
  
  ;read in redshift data
  readcol,path_to_output+'/allslits.seplines.zmosfire.qflag.txt',zmask,zslit,zobject,zredshift,format='X,A,I,I,F,X,X,X,X',/SILENT
  
  
  ;open file to write info
  openw,lun,backup_folder_name+'_info.txt',/get_lun
  printf,lun,'File created on:'
  printf,lun,systime()
  printf,lun,'left side, new... right side, old
  printf,lun
  openw,lun2,backup_folder_name+'_diff_info.txt',/get_lun
  printf,lun2,'File created on:'
  printf,lun2,systime()
  printf,lun2
  openw,lun3,backup_folder_name+'_undone_info.txt',/get_lun
  printf,lun3,'File created on:'
  printf,lun3,systime()
  printf,lun3
  
  ;open plot to plot to...
  data=readfits(filenames[0],hdr,exten_no=1,/silent)
  PS_Start, Filename=backup_folder_name+'_'+strcompress(sxpar(hdr,'MSKNM'),/remove_all)+'_comparison.ps',/quiet
  
  ;loop through 1d spectra
  for i=0,n_elements(filenames)-1 do begin
    old_mask_name=strcompress(sxpar(hdr,'MSKNM'),/remove_all)
    
    ;read in file
    temp=readfits(filenames[i],hdr0,exten_no=0,/silent)
    data=readfits(filenames[i],hdr,exten_no=1,/silent)
    ydata=readfits(filenames[i],exten_no=5,/silent)
    
    ;start new .ps file for each mask.
    if strcompress(sxpar(hdr,'MSKNM'),/remove_all) ne old_mask_name then begin
      ps_end
      PS_Start, Filename=backup_folder_name+'_'+strcompress(sxpar(hdr,'MSKNM'),/remove_all)+'_comparison.ps',/quiet
    endif
    
    ;move files to backup
    ;;;;;;;;   spawn,'mv '+filenames[i]+' '+backup_folder_name
    
    ;check AUTOEX/WBYHAND/CBYHAND/BLIND flags
    if  $ ;sxpar(hdr,'WBYHAND') eq 1 or
      sxpar(hdr,'AUTOEX')  eq 0 or $
      sxpar(hdr,'CBYHAND') eq 1 or $
      sxpar(hdr,'BLIND')   eq 1 then begin
      if not sxpar(hdr,'BLIND') then printf,lun,filenames[i]+' not done because of fits hdr flags.'
      if sxpar(hdr,'BLIND') eq 0 then printf,lun3,filenames[i],' because of fits hdr flags.'
    endif else begin
    
    
      ;check for two objects too close to each other.
    
      ;read in variables from 1d spectra
      ypos=sxpar(hdr,'YPOS')
      width=sxpar(hdr,'WIDTH')
      N_BINS=sxpar(hdr,'NBINS')
      left_bins=[]
      right_bins=[]
      FOR jj=0,N_BINS-1 do begin
        left_bins=[left_bins,sxpar(hdr,'LBINW'+ssi(jj+1))]
        right_bins=[right_bins,sxpar(hdr,'RBINW'+ssi(jj+1))]
      endfor
      wavel=(sxpar(hdr,'CRVAL1')+findgen(n_elements(data))*sxpar(hdr,'CDELT1'))
      sci_file_name=sxpar(hdr,'FILNM')
      printf,lun,filenames[i],' '+sci_file_name
      print,filenames[i]
      
      ;read in NEW! 2d image
      if file_test(path_to_output+'/'+sci_file_name) then begin
        sciimg=readfits(path_to_output+'/'+sci_file_name,D2hdr,exten_no=1,/silent)
        var_img=readfits(path_to_output+'/'+sci_file_name,D2hdr,exten_no=3,/silent)
        ;clean image
        index=where(var_img eq 99,/null)
        sciimg[index]=0.0
        var_img[index]=0.0
        index=where(finite(var_img) eq 0,/null)
        sciimg[index]=0.0
        var_img[index]=0.0
        var_img=double(var_img)
        var_img=var_img*var_img
        index=where(var_img lt 0,ct)
        var_img[INDEX]=ABS(var_img[INDEX])
        ;        index=where(var_img eq 0,ct)
        ;        if ct ne 0 then var_img[index]=9801.0
        ;print,'there are ',ct,' pixels with NEGATIVE value for  variance '
        
        
        ;convert wavel to pixels in new img.
        wavel2d=(sxpar(D2hdr,'CRVAL1')+findgen(n_elements(sciimg[*,0]))*sxpar(D2hdr,'CDELT1'))
        left_pixels=[]
        right_pixels=[]
        FOR jj=0,N_BINS-1 do begin
          index=where(abs(wavel2d-left_bins[jj]) eq min(abs(wavel2d-left_bins[jj])))
          left_pixels=[left_pixels,index[0]]
          index=where(abs(wavel2d-right_bins[jj]) eq min(abs(wavel2d-right_bins[jj])))
          right_pixels=[right_pixels,index[0]]
          printf,lun,'LBIN'+ssi(jj+1),left_pixels[jj], sxpar(hdr,'LBIN'+ssi(jj+1))
          printf,lun,'RBIN'+ssi(jj+1),right_pixels[jj],sxpar(hdr,'RBIN'+ssi(jj+1))
        endfor
        
        ;collapse columns based on wavelength
        newydata=replicate(0.0,n_elements(sciimg[0,*]))
        newydataerr=replicate(0.0,n_elements(sciimg[0,*]))
        FOR jj=0,N_BINS-1 do begin
          x0=left_pixels[jj]
          x=right_pixels[jj]
          for k=0,n_elements(newydata)-1 do begin
            if x0 lt x then index=indgen(abs(x-x0))+x0 else index=indgen(abs(x-x0))+x
            
            if sxpar(hdr,'CONTMODE') eq 1 then begin
              noiseguess=[]
              for j=0,n_elements(index)-1 do noiseguess=[noiseguess,total(var_img[index[j],*])]
              goodpix=where(noiseguess lt median(noiseguess))
              badpix=where(noiseguess ge median(noiseguess))
              index=index[goodpix]
            endif
            for j=0,n_elements(index)-1 do newydata[k]=newydata[k]+sciimg[index[j],k]
            for j=0,n_elements(index)-1 do newydataerr[k]=newydataerr[k]+var_img[index[j],k]
          endfor;k=0,n_elements(newydata)-1
        endfor;jj=0,N_BINS-1
        
        ;plot y cross sections
        ;        cgplot,ydata,title=filenames[i]+'--'+ssi(sxpar(hdr,'CONTMODE')),/nodata
        ;        cgplot,ydata,thick=6,/overplot
        
        ;MASK SPECIFIC ANOMALIES
        newydata2=newydata
        bot_slit_shift=5
        case strcompress(sxpar(hdr,'MSKNM'),/remove_all) of
          'ae1_03': begin
          
          ; if sxpar(hdr,'SLITNM') eq 1601 then ypos=ypos-bot_slit_shift
          end
          'ae1_05': begin
            delta=8
            if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'H' then delta=9
            ypos=ypos+delta
            if sxpar(hdr,'SLITNM') eq 10454 then ypos=ypos-delta
            if sxpar(hdr,'SLITNM') eq 1002 then ypos=ypos-bot_slit_shift
          end
          'ae2_03': begin
            ypos=ypos+11
            if sxpar(hdr,'SLITNM') eq 9021 then ypos=ypos-11
            if sxpar(hdr,'SLITNM') eq 1361 then ypos=ypos-bot_slit_shift
          end
          'co2_01': begin
            delta=8
            if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'K' then delta=16
            ypos=ypos+delta
            if sxpar(hdr,'SLITNM') eq 1145 then ypos=ypos-delta
            if sxpar(hdr,'SLITNM') eq 1679 then ypos=ypos-bot_slit_shift
          end
          'co2_03': begin
            delta=8
            if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'J' then delta=12
            ypos=ypos+delta
            if sxpar(hdr,'SLITNM') eq 14172 then ypos=ypos-delta
            if sxpar(hdr,'SLITNM') eq 12104 then ypos=ypos-bot_slit_shift
          end
          'co3_01': begin
          
          ;if sxpar(hdr,'SLITNM') eq 2357 then ypos=ypos-bot_slit_shift
          end
          'gn2_04': begin
            ypos=ypos+8
            if sxpar(hdr,'SLITNM') eq 4336 then ypos=ypos-8
            if sxpar(hdr,'SLITNM') eq 2182 then ypos=ypos-bot_slit_shift
          end
          'gn2_05': begin
            delta=8
            if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'H' then delta=9
            if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'K' then delta=11
            ypos=ypos+delta
            if sxpar(hdr,'SLITNM') eq 9994 then ypos=ypos-delta
            if sxpar(hdr,'SLITNM') eq 7195 then ypos=ypos-bot_slit_shift
          end
          'gs2_01': begin
            delta=5
            if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'K' then delta=9
            ypos=ypos+delta
            if sxpar(hdr,'SLITNM') eq 30874 then ypos=ypos-delta
            if sxpar(hdr,'SLITNM') eq 25572 then ypos=ypos-delta
            if sxpar(hdr,'SLITNM') eq 22627 then ypos=ypos-bot_slit_shift
          end
          'ud1_01': begin
            delta=6
            if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'H' then delta=7
            ypos=ypos+delta
            if sxpar(hdr,'SLITNM') eq 21323 then ypos=ypos-delta
            if sxpar(hdr,'SLITNM') eq 6187 then ypos=ypos-bot_slit_shift
          end
        endcase
        
        
        ;normalize cause of errors.
        newydataerr=newydataerr/(max(newydata)*max(newydata))
        newydata=newydata/max(newydata)
        
        cgplot,newydata,thick=3,color='red',title=filenames[i]+'--'+ssi(sxpar(hdr,'CONTMODE'))
        cgerrplot,findgen(n_elements(newydata)),newydata-sqrt(newydataerr),newydata+sqrt(newydataerr)
        cgplot,[ypos,ypos],[min(newydata),max(newydata)],/overplot
        
        ;        stop
        ;fit new gaussian
        xpos=ypos
        
        lower=round(xpos-mwidth)
        if lower lt 0 then lower = 0
        if lower ge n_elements(newydata) then lower = n_elements(newydata)-1
        upper=round(xpos+mwidth)
        if upper ge n_elements(newydata) then upper = n_elements(newydata)-1
        if upper lt 0 then upper = 0
        if upper ne lower then begin
          yfit=newydata[lower:upper]
          yfiterr=sqrt(newydataerr[lower:upper])
          nterms=4
          
          pi =[{fixed:0, limited:[1,1], limits:[0.D,max(newydata)*1.2]},$ ;peak value
            {fixed:0, limited:[1,1], limits:[-3.D,n_elements(newydata)+2]},$ ;peak centroid
            {fixed:0, limited:[0,0], limits:[0.D,0.D]},$ ;sigma
            {fixed:1, limited:[0,0], limits:[0.D,0.D]}];,$ ;linear bkgnd term
          estimates=[ max(yfit),xpos,2.0,0.0]
          xfit=findgen(upper-lower+1)+lower
          dummy=MPFITPEAK(double(xfit),double(yfit),$
            coeff,nterms=nterms,error=yfiterr,sigma=gauss_sigma,/gaussian,$
            estimates=estimates,parinfo=pi,status=status,chisq=chisq, ERRMSG=errmsg)
            
            
          if status eq 1 or status eq 3 then begin
            fwhm_fit=2.0*SQRT(2.0*ALOG(2.0))*coeff[2]
            if coeff[1] gt -3 and coeff[1] lt n_elements(newydata)+3 then begin
              coeff[1]=round(coeff[1])
              fwhm_fit=round(fwhm_fit)
              minw=sxpar(hdr,'MINW')
              if minw ne -1 and fwhm_fit lt minw then fwhm_fit=minw
              mom1=total(xfit*yfit)/total(yfit)
              mom2=2.355*sqrt(abs(total(xfit*xfit*yfit)/total(yfit) - $
                (total(xfit*yfit)/total(yfit))*(total(xfit*yfit)/total(yfit))))
              chisq=chisq/float(n_elements(yfit)-1.0)
              newcenter=ROUND(coeff[1])
              newwidth=ROUND(fwhm_fit)
              printf,lun,'Width: ',newwidth,width
              printf,lun,'Center: ',newcenter,ypos
              if newwidth ne width or newcenter ne ypos then begin
                printf,lun2,filenames[i]
                printf,lun2,'Width: ',newwidth,width
                printf,lun2,'Center: ',newcenter,ypos
                printf,lun2
              endif
              
              cgplot,findgen(upper-lower+1)+lower,yfit,color='green',/overplot
              cgplot,findgen(upper-lower)+lower,dummy,color='purple',/overplot
              
              if sxpar(hdr,'WBYHAND') eq 1 then newwidth = width
              
              
              if abs(newwidth - width) le 1 and abs(newcenter - ypos) le 3 then begin
                CGTEXT,0.16,0.85,'EXTRACTED',/normal
                
                ;extract and save
                bmep_extraction_simple,sciimg,var_img,newydata,newcenter,newwidth,f,ferr,fopt,fopterr,p
                
                p=p/max(p)
                p=p*max(yfit)
                
                cgplot,p,color='blue',/overplot
                
                index=where(zobject eq sxpar(hdr,'SLITNM'),ct)
                if ct eq 1 then begin
                  redshift=zredshift[index[0]]
                  wavel=wavel/(redshift+1.0)
                  wavel2d=wavel2d/(redshift+1.0)
                  count=0
                  readcol,'~/Dropbox/siana_group/bill/line_list.txt',vertlines,vertlinenames,format='F,A',/silent
                  
                  if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'K' then index=where(wavel lt 6800 and wavel gt 6500,count)
                  if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'J' then index=where(wavel lt 3900 and wavel gt 3700,count)
                  if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'H' then index=where(wavel lt 5100 and wavel gt 4800,count)
                  if count gt 10 then begin
                    cgplot,wavel[index],data[index],title=filenames[i],ytitle='red=new, black=old',$
                      yrange=[bmep_percent_cut(data[index],2),bmep_percent_cut(data[index],98)],/xs
                    for ijk=0,n_elements(vertlines)-1 do oplot_vert,vertlines[ijk],$
                      minmax(data[index]),name=vertlinenames[ijk],ynamepos=median(data[index]),charsize=1.0,/cg
                    count=0
                    if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'K' then index=where(wavel2d lt 6800 and wavel2d gt 6500,count)
                    if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'J' then index=where(wavel2d lt 3900 and wavel2d gt 3700,count)
                    if strcompress(sxpar(hdr,'FILTNM'),/remove_all) eq 'H' then index=where(wavel2d lt 5100 and wavel2d gt 4800,count)
                    
                    if count gt 10 then begin
                      cgplot,wavel2d[index],fopt[index],color='red',title=filenames[i],ytitle='red=new, black=old',$
                        yrange=[bmep_percent_cut(fopt[index],2),bmep_percent_cut(fopt[index],98)],/xs
                      for ijk=0,n_elements(vertlines)-1 do oplot_vert,vertlines[ijk],$
                        minmax(fopt[index]),name=vertlinenames[ijk],ynamepos=median(fopt[index]),charsize=1.0,/cg
                    endif;count gt 10
                  endif;count gt 10
                endif ;count gt 1
                
                
                
              ;save the damn files...
              ;                sxaddpar,hdr,'WIDTH',newwidth,/savecomment
              ;                sxaddpar,hdr,'YEXPECT',ypos,/savecomment
              ;                sxaddpar,hdr,'YPOS',newcenter,/savecomment
              ;                sxaddpar,hdr,'MOM1',mom1,/savecomment
              ;                sxaddpar,hdr,'MOM2',mom2,/savecomment
              ;                sxaddpar,hdr,'GAUSRCHI',chisq,/savecomment
              ;
              ;                sxaddpar,hdr0,'WIDTH',newwidth,/savecomment
              ;                sxaddpar,hdr0,'YEXPECT',ypos,/savecomment
              ;                sxaddpar,hdr0,'YPOS',newcenter,/savecomment
              ;                sxaddpar,hdr0,'MOM1',mom1,/savecomment
              ;                sxaddpar,hdr0,'MOM2',mom2,/savecomment
              ;                sxaddpar,hdr0,'GAUSRCHI',chisq,/savecomment
              ;
              ;                writefits,filenames[i],'',hdr0
              ;                writefits,filenames[i],fopt,hdr,/append
              ;                writefits,filenames[i],fopterr,hdr,/append
              ;                writefits,filenames[i],f,hdr,/append
              ;                writefits,filenames[i],ferr,hdr,/append
              ;
              ;                sxaddpar,hdr,'CRVAL1',1
              ;                sxaddpar,hdr,'CDELT1',1.0
              ;                sxaddpar,hdr,'CRPIX1',1
              ;                sxaddpar,hdr,'CTYPE1','LINEAR'
              ;                writefits,filenames[i],newydata,hdr,/append
                
                
              endif else begin ;new width and center OK
                printf,lun,'Bad new width or center. Too different'
                printf,lun3,filenames[i],'Bad new width or center. Too different';,newwidth,width,newcenter,ypos
                
              endelse
            endif else begin;coeff makes sense
              printf,lun,'Bad gaussian fit. (coeff out of good range.)'
              printf,lun3,filenames[i],lun,'Bad gaussian fit. (coeff out of good range.)'
            endelse
          endif else begin;fit status eq 1
            printf,lun,'Bad gaussian fit. (status ne 1)'
            printf,lun3,filenames[i],'Bad gaussian fit. (status ne 1)'
          endelse
        endif else begin;upper ne lower
          printf,lun,'Bad y position value',ypos
          printf,lun3,filenames[i],'Bad y position value',ypos
        endelse
      endif else begin ;file test
        printf,lun,'ERROR FINDING 2D FILE, '+path_to_output+'/'+sci_file_name
        printf,lun3,filenames[i],'ERROR FINDING 2D FILE, '+path_to_output+'/'+sci_file_name
      endelse
    endelse; header flags ok
    ;end loop
    printf,lun
  endfor
  
  ;close file and clean up
  close,lun
  free_lun,lun
  close,lun2
  free_lun,lun2
  close,lun3
  free_lun,lun3
  ps_end
  !p.multi=[0,1,1]
  
  bmep_mosdef_update_yexpect
  
  print,'there are ',file_lines(backup_folder_name+'_undone_info.txt')-3,' undone'
  
  print,'that took ',round(systime(/seconds)- starttime),' seconds
  cd,original_dir
  print,'done rereduciing'
end







pro bmep_mosdef_rewidth

  ;setup
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  !except=2 ;see division by zero errors instantly.
  astrolib
  starttime=systime(/seconds)
  mwidth=8 ; width to fit a gaussian to...
  !p.multi=[0,1,1]
  
  ;creat list of 1d spec
  
  if not keyword_set(path_to_output) then path_to_output='~/mosfire/output/idl_output/2D' ; no trailing slash
  cd,path_to_output+'/1d_extracted',current=original_dir
  
  ;read in starlist.txt
  readcol,path_to_output+'/1d_extracted/00_starinfo.txt',maskstar,$
  filtstar,objstar,yexpect_star,yactual_star,widthstar,sigmastar,/silent,format='A,A,A,F,F,F,F'

  ;read in list of 1d spectra
  spawn,'ls *.fc.1d.fits > tempfiles.txt' ;;for all files
;  spawn,'ls co3_01*.1d.fits > tempfiles.txt' ;for updating only one
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt'
  print,'list contains ',n_elements(filenames),' objects
  
  ;open file to write info
  openw,lun,'00_new_width_info.txt',/get_lun
;  printf,lun,'File created on:'
;  printf,lun,systime()
  printf,lun,'# Filename new_width old_width new-old new_center old_center slitloss'
;  printf,lun

  openw,lun3,'00_new_width_undone_info.txt',/get_lun
;  printf,lun3,'File created on:'
;  printf,lun3,systime()
;  printf,lun3
  
  ;open plot to plot to...
  PS_Start, Filename='00_new_width_comparison.ps',/quiet
  
  ;loop through 1d spectra
  for i=0,n_elements(filenames)-1 do begin
    
    ;read in file
    temp=readfits(filenames[i],hdr0,exten_no=0,/silent)
    data=readfits(filenames[i],hdr,exten_no=1,/silent)
    ydata=readfits(filenames[i],exten_no=5,/silent)
    
    
    ;check AUTOEX/WBYHAND/CBYHAND/BLIND flags
    if  $ ;sxpar(hdr,'WBYHAND') eq 1 or
      sxpar(hdr,'AUTOEX')  eq 0 or $
      sxpar(hdr,'CBYHAND') eq 1 or $
      sxpar(hdr,'BLIND')   eq 1 then begin
;      if ~sxpar(hdr,'BLIND') then printf,lun,filenames[i]+' not done because of fits hdr flags.'
      if sxpar(hdr,'BLIND') eq 0 then begin 
;        printf,lun3,filenames[i],' because of fits hdr flags.'
      
      
        ;SAVING THE cbyhand case.
        temp=readfits(filenames[i],hdr0,exten_no=0,/silent)
        ydata=readfits(filenames[i],hdr,exten_no=5,/silent)
        fopt=readfits(filenames[i],hdr,exten_no=1,/silent)
        fopterr=readfits(filenames[i],hdr,exten_no=2,/silent)
        f=readfits(filenames[i],hdr,exten_no=3,/silent)
        ferr=readfits(filenames[i],hdr,exten_no=4,/silent)
        
        sxaddpar, hdr, 'SLITLOSS', bmep_calc_slitloss(sxpar(hdr,'WIDTH')/2.355), ' Flux in slit / Total flux assuming 2D Gaussians'
        sxaddpar, hdr, 'GWIDTH', -99, ' Sigma of gaussian fit'
        sxaddpar, hdr, 'GCENT', -99, ' Central pixel of gaussian fit'
        sxaddpar, hdr, 'GAMP', -99, ' Amplitude of gaussian fit'
        sxaddpar, hdr, 'GLINEAR', -99, ' Linear continuum of gaussian fit'
        sxaddpar, hdr, 'GWIDTH_E', -99, ' 1sigma Error in Sigma of gaussian fit'
        sxaddpar, hdr, 'GCENT_E', -99, ' 1sigma Error in CENTRAL pixel of gaussian fit'
        sxaddpar, hdr, 'GAMP_E', -99,  ' 1sigma Error in Amplitude of gaussian fit'
        sxaddpar, hdr, 'GLINE_E', -99, ' 1sigma Error in Linear continuum of gaussian fit' 
        
;        filenames[i]='rw__'+filenames[i]
;        writefits,savename,'',hdr0
;        writefits,savename,fopt,hdr,/append
;        writefits,savename,fopterr,hdr,/append
;        writefits,savename,f,hdr,/append
;        writefits,savename,ferr,hdr,/append
;        writefits,savename,ydata,hdr,/append
        endif ; is not blind.
      
      
    endif else begin
    
    
      ;check for two objects too close to each other.
    
      ;read in variables from 1d spectra
      ypos=sxpar(hdr,'YPOS')
      width=sxpar(hdr,'WIDTH')
      N_BINS=sxpar(hdr,'NBINS')
      left_bins=[]
      right_bins=[]
      FOR jj=0,N_BINS-1 do begin
        left_bins=[left_bins,sxpar(hdr,'LBINW'+ssi(jj+1))]
        right_bins=[right_bins,sxpar(hdr,'RBINW'+ssi(jj+1))]
      endfor
      wavel=(sxpar(hdr,'CRVAL1')+findgen(n_elements(data))*sxpar(hdr,'CDELT1'))
      sci_file_name=sxpar(hdr,'FILNM')
;      printf,lun,filenames[i],' '+sci_file_name
      print,filenames[i]
      
      ;read in NEW! 2d image
      if file_test(path_to_output+'/'+sci_file_name) then begin
        sciimg=readfits(path_to_output+'/'+sci_file_name,D2hdr,exten_no=1,/silent)
        var_img=readfits(path_to_output+'/'+sci_file_name,D2hdr,exten_no=3,/silent)
        ;clean image
        index=where(var_img eq 99,/null)
        sciimg[index]=0.0
        var_img[index]=0.0
        index=where(finite(var_img) eq 0,/null)
        sciimg[index]=0.0
        var_img[index]=0.0
        var_img=double(var_img)
        var_img=var_img*var_img
        index=where(var_img lt 0,ct)
        var_img[INDEX]=ABS(var_img[INDEX])
        ;        index=where(var_img eq 0,ct)
        ;        if ct ne 0 then var_img[index]=9801.0
        ;print,'there are ',ct,' pixels with NEGATIVE value for  variance '
        
        
        ;convert wavel to pixels in new img.
        wavel2d=(sxpar(D2hdr,'CRVAL1')+findgen(n_elements(sciimg[*,0]))*sxpar(D2hdr,'CDELT1'))
        left_pixels=[]
        right_pixels=[]
        FOR jj=0,N_BINS-1 do begin
          index=where(abs(wavel2d-left_bins[jj]) eq min(abs(wavel2d-left_bins[jj])))
          left_pixels=[left_pixels,index[0]]
          index=where(abs(wavel2d-right_bins[jj]) eq min(abs(wavel2d-right_bins[jj])))
          right_pixels=[right_pixels,index[0]]
;          printf,lun,'LBIN'+ssi(jj+1),left_pixels[jj], sxpar(hdr,'LBIN'+ssi(jj+1))
;          printf,lun,'RBIN'+ssi(jj+1),right_pixels[jj],sxpar(hdr,'RBIN'+ssi(jj+1))
        endfor
        
        ;collapse columns based on wavelength
        newydata=replicate(0.0,n_elements(sciimg[0,*]))
        newydataerr=replicate(0.0,n_elements(sciimg[0,*]))
        FOR jj=0,N_BINS-1 do begin
          x0=left_pixels[jj]
          x=right_pixels[jj]
          for k=0,n_elements(newydata)-1 do begin
            if x0 lt x then index=indgen(abs(x-x0))+x0 else index=indgen(abs(x-x0))+x
            
            if sxpar(hdr,'CONTMODE') eq 1 then begin
              noiseguess=[]
              for j=0,n_elements(index)-1 do noiseguess=[noiseguess,total(var_img[index[j],*])]
              goodpix=where(noiseguess lt median(noiseguess))
              badpix=where(noiseguess ge median(noiseguess))
              index=index[goodpix]
            endif
            for j=0,n_elements(index)-1 do newydata[k]=newydata[k]+sciimg[index[j],k]
            for j=0,n_elements(index)-1 do newydataerr[k]=newydataerr[k]+var_img[index[j],k]
          endfor;k=0,n_elements(newydata)-1
        endfor;jj=0,N_BINS-1
        
;        newydata=ydata

        ;normalize cause of errors.
        newydataerr=newydataerr/(max(newydata)*max(newydata))
        newydata=newydata/max(newydata)
        

        
        ;        stop
        ;fit new gaussian
        xpos=ypos
        
        lower=round(xpos-mwidth)
        if lower lt 0 then lower = 0
        if lower ge n_elements(newydata) then lower = n_elements(newydata)-1
        upper=round(xpos+mwidth)
        if upper ge n_elements(newydata) then upper = n_elements(newydata)-1
        if upper lt 0 then upper = 0
        if upper ne lower then begin
          yfit=newydata[lower:upper]
          yfiterr=sqrt(newydataerr[lower:upper])
          nterms=4
          
          pi =[{fixed:0, limited:[1,1], limits:[0.D,max(newydata)*1.2]},$ ;peak value
            {fixed:0, limited:[1,1], limits:[-3.D,n_elements(newydata)+2]},$ ;peak centroid
            {fixed:0, limited:[0,0], limits:[0.D,0.D]},$ ;sigma
            {fixed:1, limited:[0,0], limits:[0.D,0.D]}];,$ ;linear bkgnd term
          estimates=[ max(yfit),xpos,2.0,0.0]
          xfit=findgen(upper-lower+1)+lower
          dummy=MPFITPEAK(double(xfit),double(yfit),$
            coeff,nterms=nterms,error=yfiterr,sigma=gauss_sigma,/gaussian,$
            estimates=estimates,parinfo=pi,status=status,chisq=chisq, ERRMSG=errmsg)
            
         chisq=chisq/float(n_elements(yfit)-1.0)
            
          if status eq 1 or status eq 3 then begin
            fwhm_fit=2.0*SQRT(2.0*ALOG(2.0))*coeff[2]
            if coeff[1] gt -3 and coeff[1] lt n_elements(newydata)+3 then begin
              
;              minw=sxpar(hdr,'MINW')
              index=where(SSS(maskstar) eq sss(sxpar(hdr,'MSKNM')) and sss(filtstar) eq sss(sxpar(hdr,'FILTNM')),ct)
              if ct ne 0 then begin
                minw=min(widthstar[index],sub)
              endif else message,'no star failure'
  

              if (minw gt 0 and fwhm_fit lt minw)  then begin 
                fwhm_fit=minw
                poststr='minw=FIXED'
                endif else poststr='minw=OK'
              
              IF sxpar(hdr,'ISSTAR') EQ 1 THEN BEGIN
                index=where(SSS(maskstar) eq sss(sxpar(hdr,'MSKNM')) and sss(filtstar) eq sss(sxpar(hdr,'FILTNM')) $
                  AND sss(objstar) eq sss(sxpar(hdr,'SLITNM')),ct)
                if ct eq 1 then fwhm_fit=widthstar[index[0]] else message,'wat'
                endif
              

              chisq=chisq/float(n_elements(yfit)-1.0)
              newcenter=coeff[1]
              newwidth=fwhm_fit

;              stop
              if sxpar(hdr,'WBYHAND') eq 1 then newwidth = width
              if sxpar(hdr,'WBYHAND') eq 1 and abs(newcenter - ypos) gt 1.0 then newcenter = ypos
;              if sxpar(hdr,'ISSTAR') then begin
              
              printf,lun,filenames[i],newwidth,width,$
                newwidth-width,newcenter,ypos,$
                bmep_calc_slitloss(fwhm_fit/(2.0*SQRT(2.0*ALOG(2.0)))),format='(A30,6F12.5)'
              
;                CGTEXT,0.56,0.86,'wbyhand: '+ssi(sxpar(hdr,'WBYHAND')),/normal
;                 CGTEXT,0.56,0.76,poststr,/normal
              if abs(newwidth - width) le 1.0 and abs(newcenter - ypos) le 1.0 then begin
;                CGTEXT,0.16,0.86,'NEW: '+ssf(newwidth)+' OLD:'+ssf(width),/normal
;                CGTEXT,0.16,0.76,'NEW:  '+ssf(newcenter)+' OLD: '+ssf(YPOS),/normal
                
                ;make a residual plot
                ;double(xfit),double(yfit)
                
                  !P.multi=[0,1,1]
                  cgplot,100*newydata/total(dummy),thick=3,color='red',$
                    title=filenames[i],xr=[ypos-width*3.0,ypos+width*3.0];+'--'+ssi(sxpar(hdr,'CONTMODE'))
                  cgerrplot,findgen(n_elements(newydata)),newydata-sqrt(newydataerr),newydata+sqrt(newydataerr)
                  cgplot,[ypos,ypos],[min(newydata),max(newydata)],/overplot
                  cgplot,findgen(upper-lower+1)+lower,100*yfit/total(dummy),color='green',/overplot
                  cgplot,xfit,100*dummy/total(dummy),color='blue',/overplot
                  cgplot,xfit,100*dummy/total(dummy)-100*yfit/total(dummy),title='fit - data',xr=[ypos-width*3.0,ypos+width*3.0]
                  !P.multi=[0,1,1]
;                  endif
                
                
                
                ;SAVING THE NORMAL CASE
                ;read in ex=0 hdr
                temp=readfits(filenames[i],hdr0,exten_no=0,/silent)
;                ydata=readfits(filenames[i],hdr,exten_no=5,/silent)
;                fopt=readfits(filenames[i],hdr,exten_no=1,/silent)
;                fopterr=readfits(filenames[i],hdr,exten_no=2,/silent)
;                f=readfits(filenames[i],hdr,exten_no=3,/silent)
;                ferr=readfits(filenames[i],hdr,exten_no=4,/silent)
                
                ;reextract
;                bmep_extraction_simple,sciimg,var_img,newydata,round(newcenter),round(newwidth),$
;                  f,ferr,fopt,fopterr,p
                
                ;add new parameters
                sxaddpar, hdr, 'SLITLOSS', bmep_calc_slitloss(newwidth/2.0*SQRT(2.0*ALOG(2.0))), $
                  ' Flux in slit / Total flux assuming 2D Gaussians'
                sxaddpar, hdr, 'GWIDTH', coeff[2], ' Sigma of gaussian fit'
                sxaddpar, hdr, 'GCENT', coeff[1], ' Central pixel of gaussian fit'
                sxaddpar, hdr, 'GAMP', coeff[0], ' Amplitude of gaussian fit'
                sxaddpar, hdr, 'GLINEAR', coeff[3], ' Linear continuum of gaussian fit'
                sxaddpar, hdr, 'GWIDTH_E', gauss_sigma[2], ' 1sigma Error in Sigma of gaussian fit'
                sxaddpar, hdr, 'GCENT_E', gauss_sigma[1], ' 1sigma Error in CENTRAL pixel of gaussian fit'
                sxaddpar, hdr, 'GAMP_E', gauss_sigma[0],  ' 1sigma Error in Amplitude of gaussian fit'
                sxaddpar, hdr, 'GLINE_E', gauss_sigma[3], ' 1sigma Error in Linear continuum of gaussian fit' 
                sxaddpar, hdr, 'YPOS', round(newcenter),' Y position of center of extraction'
                sxaddpar, hdr, 'WIDTH', round(newwidth),' Half width in n_pixels'

                
;                if ~sxpar(hdr,'ISSTAR') then begin
;                  filenames[i]='rw__'+filenames[i]
;;                  writefits,savename,'',hdr0
;;                  writefits,savename,fopt,hdr,/append
;;                  writefits,savename,fopterr,hdr,/append
;;                  writefits,savename,f,hdr,/append
;;                  writefits,savename,ferr,hdr,/append
;;                  writefits,savename,ydata,hdr,/append
;;                  writefits,savename,newydataerr,hdr,/append
;                  endif

              endif else begin ;new width and center OK
;                printf,lun,filenames[i],' Bad new width or center. Too different'
                printf,lun3,filenames[i],' Bad new width or center. Too different';,newwidth,width,newcenter,ypos
;                CGTEXT,0.16,0.96,'BAD',/normal
;                CGTEXT,0.16,0.86,'NEW: '+ssf(newwidth)+' OLD:'+ssf(width),/normal
;                CGTEXT,0.16,0.76,'NEW:  '+ssf(newcenter)+' OLD: '+ssf(YPOS),/normal
              endelse
            endif else begin;coeff makes sense
;              printf,lun,filenames[i],' Bad gaussian fit. (coeff out of good range.)'
              printf,lun3,filenames[i],lun,' Bad gaussian fit. (coeff out of good range.)'
            endelse
          endif else begin;fit status eq 1
          
;            printf,lun,filenames[i],' Bad gaussian fit. (status ne 1)'
            printf,lun3,filenames[i],' Bad gaussian fit. (status ne 1)'
          endelse
        endif else begin;upper ne lower
;          printf,lun,filenames[i],' Bad y position value',ypos
          printf,lun3,filenames[i],' Bad y position value',ypos
        endelse
      endif else begin ;file test
        printf,lun,filenames[i],' ERROR FINDING 2D FILE, '+path_to_output+'/'+sci_file_name
        printf,lun3,filenames[i],' ERROR FINDING 2D FILE, '+path_to_output+'/'+sci_file_name
      endelse
    endelse; header flags ok
    ;end loop
;    printf,lun
  endfor
  
  ;close file and clean up
  close,lun
  free_lun,lun
  close,lun3
  free_lun,lun3
  ps_end
  !p.multi=[0,1,1]
  
;  bmep_mosdef_update_yexpect
  
  print,'that took ',round(systime(/seconds)- starttime),' seconds
  cd,original_dir
  print,'done rewidthing'
end








;6544 7428m
pro bmep_mosdef_getinfo,doblind=doblind
  !except=2 ;see division by zero errors instantly.
  astrolib
  
  if not keyword_set(path_to_output) then path_to_output=$
    '~/mosfire/output/idl_output/2D/1d_extracted' ; no trailing slash
  cd,path_to_output,current=original_dir
  
  spawn,'ls *.1d.fits > filenames.txt'
  readcol,'filenames.txt',filenames,format='A'
  spawn,'rm filenames.txt'
  maskarr=[]
  filtarr=[]
  slitarr=[]
  widtharr=[]
  yposarr=[]
  yexpectarr=[]
  objnumarr=[]
  isstararr=[]
  wbyhandarr=[]
  yexpectarr=[]
  minwarr=[]
  blindarr=[]
  
  for i=0,n_elements(filenames)-1 do begin
    data=readfits(filenames[i],hdr,exten_no=1,/silent)
    substrings=strsplit(filenames[i],'.',/extract)
    if substrings[0] ne 'blind' then begin
      maskarr=[maskarr,sxpar(hdr,'MSKNM')]
      filtarr=[filtarr,sxpar(hdr,'FILTNM')]
      slitarr=[slitarr,sxpar(hdr,'SLITNM')]
      widtharr=[widtharr,fix(sxpar(hdr,'WIDTH'))]
      yposarr=[yposarr,sxpar(hdr,'YPOS')]
      objnumarr=[objnumarr,sxpar(hdr,'OBJNUM')]
      isstararr=[isstararr,sxpar(hdr,'ISSTAR')]
      wbyhandarr=[wbyhandarr,sxpar(hdr,'WBYHAND')]
      yexpectarr=[yexpectarr,sxpar(hdr,'YEXPECT')]
      minwarr=[minwarr,sxpar(hdr,'MINW')]
      blindarr=[blindarr,sxpar(hdr,'BLIND')]
    endif
  endfor
  
  
  index=bsort(filtarr)
  maskarr=maskarr[index]
  filtarr=filtarr[index]
  slitarr=slitarr[index]
  widtharr=widtharr[index]
  yposarr=yposarr[index]
  objnumarr=objnumarr[index]
  isstararr=isstararr[index]
  wbyhandarr=wbyhandarr[index]
  yexpectarr=yexpectarr[index]
  minwarr=minwarr[index]
  blindarr=blindarr[index]
  
  index=bsort(objnumarr)
  maskarr=maskarr[index]
  filtarr=filtarr[index]
  slitarr=slitarr[index]
  widtharr=widtharr[index]
  yposarr=yposarr[index]
  objnumarr=objnumarr[index]
  isstararr=isstararr[index]
  wbyhandarr=wbyhandarr[index]
  yexpectarr=yexpectarr[index]
  minwarr=minwarr[index]
  blindarr=blindarr[index]
  
  index=bsort(slitarr)
  maskarr=maskarr[index]
  filtarr=filtarr[index]
  slitarr=slitarr[index]
  widtharr=widtharr[index]
  yposarr=yposarr[index]
  objnumarr=objnumarr[index]
  isstararr=isstararr[index]
  wbyhandarr=wbyhandarr[index]
  yexpectarr=yexpectarr[index]
  minwarr=minwarr[index]
  blindarr=blindarr[index]
  
  ;    index=sort(filtarr)
  ;  maskarr=maskarr[index]
  ;  filtarr=filtarr[index]
  ;  slitarr=slitarr[index]
  ;  widtharr=widtharr[index]
  ;  yposarr=yposarr[index]
  ;  objnumarr=objnumarr[index]
  ;  isstararr=isstararr[index]
  ;  wbyhandarr=wbyhandarr[index]
  ;  yexpectarr=yexpectarr[index]
  
  index=bsort(maskarr)
  maskarr=maskarr[index]
  filtarr=filtarr[index]
  slitarr=slitarr[index]
  widtharr=widtharr[index]
  yposarr=yposarr[index]
  objnumarr=objnumarr[index]
  isstararr=isstararr[index]
  wbyhandarr=wbyhandarr[index]
  yexpectarr=yexpectarr[index]
  minwarr=minwarr[index]
  blindarr=blindarr[index]
  
  
  
  
  openw,lun,'00_extract_info.txt',/get_lun
  printf,lun,'mask filter slit-#[*] width[*](starwidth) ypos yposdifference[!!!]'
  printf,lun,'* after object number means this is the STAR'
  printf,lun,'* after width means the width was edited by HAND'
  printf,lun,'& after width means blindly extracted'
  printf,lun,'!!! after yposdiff means this value is more than 4 different than others'
  print,     'mask filter slit-#[*] width[*](starwidth) ypos yposdifference'
  print,     '* after object number means this is the STAR'
  print,     '* after width means the width was edited by HAND'
  print,     '& after width means blindly extracted'
  print,     '!!! after yposdiff means this value is more than 4 different than others'
  for i=0,n_elements(slitarr)-1 do begin
  
    ;space out different objects
    if i gt 0 and slitarr[i] ne slitarr[i-1] then print
    if i gt 0 and slitarr[i] ne slitarr[i-1] then printf,lun
    
    ;space out masks
    if i gt 0 and maskarr[i] ne maskarr[i-1] then begin
      ;      printf,lun
      ;      print
      printf,lun,'======================================='
      print,'======================================='
      printf,lun
      print
    endif
    
    if isstararr[i]  eq 1 then suffix='*'  else suffix=' '
    if wbyhandarr[i] eq 1 then suffix2='*' else $
      if blindarr[i]   eq 1 then suffix2='&' else suffix2=' '
    if (i gt 0) and $
      (slitarr[i] eq slitarr[i-1]) and $
      (objnumarr[i] eq objnumarr[i-1]) and $
      (isstararr[i] ne 1) and $
      (abs((yposarr[i]-yexpectarr[i]) - (yposarr[i-1]-yexpectarr[i-1])) gt 4.0) $
      then suffix3='!!!' else suffix3='   '
      
    print,     maskarr[i],filtarr[i],slitarr[i],objnumarr[i],$
      widtharr[i],minwarr[I],yposarr[i],yposarr[i]-yexpectarr[i],$
      format='(A7,A2,I8,"-",I1,"'+suffix+'",I3,"'+suffix2+'(",I1,") ",I3,F7.1,"'+suffix3+'")
    printf,lun,maskarr[i],filtarr[i],slitarr[i],objnumarr[i],$
      widtharr[i],minwarr[I],yposarr[i],yposarr[i]-yexpectarr[i],$
      format='(A7,A2,I8,"-",I1,"'+suffix+'",I3,"'+suffix2+'(",I1,") ",I3,F7.1,"'+suffix3+'")
  ;  stop
  endfor
  index=where(wbyhandarr eq 1)
  print,n_elements(index)
  
  close,lun
  free_lun,lun
  
  cd,original_dir
end

pro bmep_mosdef_blind_extract,yexpect,width,ny,ptemp,sciimg,var_img, $
    $ ;OUTPUTS
    f,ferr,fopt,fopterr,p
  f=[]
  ferr=[]
  fopt=[]
  fopterr=[]
  
  bottomint=fix(yexpect-width)
  if bottomint lt 0 then bottomint=0
  topint=fix(yexpect+width)
  if topint gt ny-1 then topint=ny-1
  
  ;shift p using magic...
  xarr_small=findgen(topint-bottomint+1)+bottomint ;create new xarr_small
  p=replicate(0.0,n_elements(sciimg[0,*]))
  bpm=replicate(1.0,n_elements(sciimg[0,*]))
  p[xarr_small]=ptemp
  if n_elements(xarr_small) ne n_elements(ptemp) then print,'ERROR: WRONG NUMBER OF ELEMENTS SOMEWHERE'
  
  for j=0,n_elements(sciimg[*,0])-1 do begin
    f=[f,total(sciimg[j,bottomint:topint])]
    ferr=[ferr,sqrt(total(var_img[j,bottomint:topint]))]
    
    bpm=replicate(1.0,n_elements(p))
    
    var_opt=var_img[j,*]
    index=where(var_opt[xarr_small] eq 0.0 or var_opt[xarr_small] eq 99.00,/null)
    bpm[xarr_small[index]]=0.0
    var_opt[xarr_small[index]]=1.0
    
    ;calc optimal flux as in horne
    numerator=total(bpm[xarr_small]*p[xarr_small]*$
      (sciimg[j,xarr_small])/var_opt[xarr_small])
    denominator=total(bpm[xarr_small]*p[xarr_small]*p[xarr_small]/var_opt[xarr_small])
    if denominator ne 0 then flux_opt=numerator/denominator else flux_opt=0.0
    
    ;calc optimal flux error as in horne.
    numerator=total(bpm[xarr_small]*p[xarr_small]) ;
    denominator=total(bpm[xarr_small]*p[xarr_small]*p[xarr_small]/var_opt[xarr_small])
    if denominator ne 0 then err_opt=sqrt(abs(numerator/denominator)) else err_opt=0.0
    
    
    fopt=[fopt,flux_opt]
    fopterr=[fopterr,err_opt]
  endfor ; cols of data! - end of the optimal section of extraction.
  
  
end

pro bmep_mosdef_blind_save,savepath,maskname,filtername,slitname,$
    objnum,extrainfo1,extrainfo2,yexpect,width,isstar,$
    f,ferr,fopt,fopterr,p,min_width
  ;  prefix='blind.'
  prefix=''
  suffix=''
  ;  print,'object number, ',objnum
  doSave=0
  if objnum ne 1 then suffix='.'+ssi(objnum)

  thefilename=savepath+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits'
  
  if file_test(thefilename) eq 0 then doSave=1 $
  else begin
    data=readfits(thefilename,hdr,exten_no=1,/silent)
    IF sxpar(hdr,'BLIND') EQ 1 THEN doSave=1 else doSave=0
    endelse

if doSave eq 1 then begin
  ;    print,'WOULD HAVE SAVED: '+prefix+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits'
  header=bmep_blind_hdr('',extrainfo1,extrainfo2,yexpect,width,isstar,objnum,min_width,/exten)
  writefits,thefilename,'',header
  header=bmep_blind_hdr(fopt,extrainfo1,extrainfo2,yexpect,width,isstar,objnum,min_width,/image)
  writefits,thefilename,fopt,header,/append
  header=bmep_blind_hdr(fopterr,extrainfo1,extrainfo2,yexpect,width,isstar,objnum,min_width,/image)
  writefits,thefilename,fopterr,header,/append
  header=bmep_blind_hdr(f,extrainfo1,extrainfo2,yexpect,width,isstar,objnum,min_width,/image)
  writefits,thefilename,f,header,/append
  header=bmep_blind_hdr(ferr,extrainfo1,extrainfo2,yexpect,width,isstar,objnum,min_width,/image)
  writefits,thefilename,ferr,header,/append
  header=bmep_blind_hdr(p,extrainfo1,extrainfo2,yexpect,width,isstar,objnum,min_width,/image,/no_wave_info)
  writefits,thefilename,p,header,/append
endif

;print,maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits'


;if doSave eq 0 then print,'not saved ' else print,'saved'

end



pro bmep_mosdef_blind,path_to_dropbox=path_to_dropbox,path_to_output=path_to_output
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  starttime=systime(/seconds)
  norepeat=0 ; 1 skips objects already done, 0 redoes objects done blindly.
  !except=2 ;see division by zero errors instantly.
  astrolib
  !p.multi=[0,2,2]
  width=5
  
  if not keyword_set(path_to_output) then path_to_output='~/mosfire/output/idl_output/2D' ; no trailing slash
  cd,path_to_output,current=original_dir
  
  ;parse folder into different masks!!
  spawn,'ls *.2d.fits > tempfiles.txt'
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt'
  
  ;parse names
  masks=[]
  filters=[]
  objects=[]
  
  for i=0,n_elements(filenames)-1 do begin
    substrings=strsplit(filenames[i],'.',/extract)
    if n_elements(substrings) eq 5 then begin
      masks=[masks,substrings[0]]
      filters=[filters,substrings[1]]
      objects=[objects,substrings[2]]
    endif
  end
  
  savepath=path_to_output+'/1d_extracted/'
  
  ;find 1d extractions of non-primary objects
  print,'creating secondary obj info'
  spawn,'ls '+savepath+'*.1d.fits > tempfiles.txt'
  readcol,'tempfiles.txt',filenames1d,format='A',/silent
  spawn,'rm tempfiles.txt'
  
  openw,lun,savepath+'00_be_info.txt',/get_lun
  for i=0,n_elements(filenames1d)-1 do begin
    data=readfits(filenames1d[i],hdr,exten_no=1,/silent)
    objnum=sxpar(hdr,'OBJNUM')
    if objnum ne 1 and sxpar(hdr,'BLIND') ne 1 then begin
      mask=sxpar(hdr,'MSKNM')
      filter=sxpar(hdr,'FILTNM')
      slit=sxpar(hdr,'SLITNM')
      width=FLOAT(round(sxpar(hdr,'WIDTH')))
      width_sigma=sxpar(hdr,'GWIDTH')
      minw=float(sxpar(hdr,'MINW'))
      if width_sigma lt minw then width_sigma=minw
      ypos=sxpar(hdr,'YPOS')
      yexpect=sxpar(hdr,'YEXPECT')
      w_actual_squared=width_sigma*width_sigma-minw*minw
      if wscale lt 1.0 then wscale=1.0
      printf,lun,mask,filter,slit,objnum,width,$
        w_actual_squared,ypos,yexpect, ypos-yexpect,format='(A10,A4,I10,I3,I3,F12.5,I4,F9.2,F9.2)'
    endif
  endfor;n_ele filenames1d
  close,lun
  free_lun,lun
  print,'done creating secondary obj info'
  
  
  
  
  readcol,savepath+'00_be_info.txt',npmaskarr,npfilterarr,npslitarr,npobjnumarr,$
    npwidtharr,w_actual_squared,npyposarr,npyexpectarr,npyshiftarr,$
    format="A,A,I,I,I,F,I,F,F"
    
    
  if norepeat eq 0 then PS_Start, Filename=savepath+'00_blind_comparison.ps',/quiet
  
  for i=0,n_elements(objects)-1 do begin
    slitname=objects[i]
    maskname=masks[i]
    filtername=filters[i]
    filename=maskname+'.'+filtername+'.'+slitname+'.2d.fits'
    
    ;  print,filename
    
    ;read in files
    sciimg=readfits(filename,shdr, /SILENT,exten_no=1)
    sciimg=double(sciimg)
    
    index=where(finite(sciimg) eq 0,/null)
    sciimg[index]=0.0
    
    ;calculate variance image
    var_img=readfits(filename, /SILENT,exten_no=3)
    
    ;clean image
    index=where(var_img eq 99,/null)
    sciimg[index]=0.0
    var_img[index]=0.0
    index=where(finite(var_img) eq 0,/null)
    sciimg[index]=0.0
    var_img[index]=0.0
    index=where(var_img lt 0,ct,/null)
    var_img=double(var_img)
    var_img=var_img*var_img ;actually was sigma image
    
    
    ;calculate wavel.
    wavel=(sxpar(shdr,'CRVAL1')+findgen(n_elements(sciimg[*,0]))*sxpar(shdr,'CDELT1'))
    ny=n_elements(sciimg[0,*])
    nx=n_elements(sciimg[*,0])
    
    ;calculate where object SHOULD be!
    slitlistfile=path_to_output+'/00mask_info/'+maskname+'_SlitList.txt'
    if ~file_test(slitlistfile) then message,'SLITLIST FILE MISSING... '+SLITLISTFILE
    
    ;1 2 17  45.61 -5  10  9.37  0.70  7.01  23763 800.00  -0.38 2 17  45.59 -5  10  9.65
    readcol,slitlistfile,slitnamearr,priorityarr,offsetarr,format='X,X,X,X,X,X,X,X,X,I,F,F,X,X,X,X,X,X',/silent
    index=where(slitnamearr eq slitname,ct)
    yexpect=-1
    isstar=-1
    yshift=0.0
    if ct eq 1 then begin
      pixscale=0.1799
      midpoint=ny/2
      yexpect=midpoint+offsetarr[index[0]]/pixscale
      ;check if object is a star
      if abs(priorityarr[index]) eq 1 then isstar=1
      ;        isstar=0
      readcol,path_to_output+'/1d_extracted/00_starinfo.txt',maskstar,$
        filtstar,objstar,yexpect_star,yactual_star,widthstar,sigmastar,/silent,format='A,A,A,F,F,F,F'
      index=where(maskstar eq maskname and filtstar eq filtername,ct)
      if ct ne 0 then begin
        ;          print,ct,' number of stars found for ',maskname,' ',filtername
        yshift=avg(yexpect_star[index] - yactual_star[index])
        yexpect=yexpect-yshift
        yexpect=round(yexpect)
        width=fix(min(widthstar[index],sub))
        starwidth=width
        min_width=starwidth
        yposstar=yactual_star[index[sub]]
        starfile=maskname+'.'+filtername+'.'+objstar[index[sub]]+'.1d.fits'
        p=readfits('1d_extracted/'+starfile,exten_no=5,/silent)
        
        ;lower extraction limit
        bottomint=fix(yposstar-width)
        if bottomint lt 0 then bottomint=0
        
        ;upper extraction limit
        topint=fix(yposstar+width)
        if topint gt n_elements(p)-1 then topint=n_elements(p)-1
        
        ;calculate array to extract
        xarr_small=findgen(topint-bottomint+1)+bottomint
        xarr_big  =findgen(n_elements(p))
        
        index=where(xarr_big lt bottomint or xarr_big gt topint,ct)
        if ct ne 0 then p[index]=0.0
        
        gresult=MPFITPEAK(double(xarr_big),double(p),$
          coeff,nterms=3,/gaussian)
        p=gresult
        
        ;set p to 0 outside the width of extraction.
        index=where(xarr_big lt bottomint or xarr_big gt topint,/null)
        p[index]=0.0
        
        ;normalize p
        index=where(p lt 0,/null)
        p[index]=0.0
        if total(p) eq 0 then print,'ERROR: P IS 0 EVERYWHERE'
        p=p/total(p)
        
        ;create temp p of just the gaussian.
        ptemp=p[xarr_small] ; uses OLD xarr_small
        
        
      endif else begin
        print,'no matching stars found ',maskname,' ',filtername
        yexpect=-1
      endelse
      
      
      
      
      ;IF THE OBJECT IS EXTRACTED IN OTHER BANDS, THEN USE THEIR WIDTH, NOT THE STAR.
      index=where(npmaskarr eq maskname and npslitarr eq slitname AND npobjnumarr eq 1,ct)
      if ct gt 0 then begin
        
        endif;fixing objects 
      
      
      
      
      
      
      
      
      
      
      ;      print,'slitname, yexpect, midpoint, yshift, width'
      print,maskname,' ', filtername,' ', slitname,' ',1, yexpect, midpoint, yshift, width
      
    endif else print,'no object found in the slitlist file?!?!?!?'
    
    
    ;calculate info to add to hdr
    ;its a 2xn array where n is number of things
    extrainfo1=[$
      'CRVAL1',$
      'CDELT1',$
      'CRPIX1',$
      'CTYPE1',$
      'EXPTIME',$
      'FILNM',$
      'MSKNM',$
      'FILTNM',$
      'SLITNM',$
      'ISSTAR',$
      'YEXPECT'$
      ]
      
    extrainfo2=[$
      string(sxpar(shdr,'CRVAL1')),$
      string(sxpar(shdr,'CDELT1')),$
      string(sxpar(shdr,'CRPIX1')),$
      'LINEAR',$
      string(sxpar(shdr,'EXPTIME')),$
      filename,$
      maskname,$
      filtername,$
      slitname, $
      ssi(isstar), $
      ssf(yexpect) $
      ]
    ;comments
    extrainfo3=[$
      ' ',$
      ' ',$
      ' ',$
      ' ',$
      'total exposure time (seconds)',$
      'name of file',$
      'name of mask',$
      'name of filter',$
      'name of slit', $
      'Flag if is a star', $
      'expected y position' $
      ]
      
      
      
      
      
      
      
    yexpect=float(yexpect)
    if yexpect ne -1 then begin
    
      if slitname eq 6187 then yexpect=yexpect-5
      if slitname eq 1002 then yexpect=yexpect-5
      
      
      bmep_mosdef_blind_extract,yexpect,width,ny,ptemp,sciimg,var_img, $
        $ ;OUTPUTS
        f,ferr,fopt,fopterr,p
      objnum=1
      bmep_mosdef_blind_save,savepath,maskname,filtername,slitname,$
        objnum,extrainfo1,extrainfo2,yexpect,width,isstar,$
        f,ferr,fopt,fopterr,p,min_width
        
      ;plot a comparison if needed...
      if norepeat eq 0 then begin
        suffix='' ;ssi(objnum)
        ;search for existing 1d file
        if file_test(savepath+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits')  then begin
          data=readfits(savepath+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits',shdr,exten_no=1,/silent)
          wavel=(sxpar(shdr,'CRVAL1')+findgen(n_elements(sciimg[*,0]))*sxpar(shdr,'CDELT1'))
          data=readfits(savepath+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits',shdr,exten_no=5,/silent)
          cgplot,data,title=maskname+'.'+filtername+'.'+slitname+suffix+' y profile',ytitle='P'
          cgplot,(p/max(p))*max(data),color='red',/overplot
          
        ;          !p.multi=[0,1,1]
        endif;file test
      endif ;norepeat eq 0 /file test
      
      ;search for np objects
      ;      ,npmaskarr,npfilterarr,npslitarr,npobjnumarr,$
      ;      npwidtharr,npwscalearr,npyposarr,npyexpectarr,npyshiftarr,
      index=where(npmaskarr eq maskname and npslitarr eq slitname,ct)
      if ct gt 0 then begin
        ;create array of objects found.
        npobjnums=rem_dup(npobjnumarr[index])
        for k=2,6 do begin
          index2=where(npobjnumarr[index] eq k,ct)
          if ct ne 0 then begin ;message,'insanity. probably an object #3 is there with no object #2'
            objnum=k
            npyexpect=round(yexpect+avg(npyshiftarr[index[index2]]))
            npwidth=round(starwidth*avg(npwscalearr[index[index2]]))
            
            print,maskname,' ', filtername,' ', slitname,' ',k, npyexpect, midpoint, yshift, npwidth
            
            ;calculate new ptemp. needed because widths are different.
            bottomint=fix(npyexpect-npwidth)
            if bottomint lt 0 then bottomint=0
            topint=fix(npyexpect+npwidth)
            if topint gt ny-1 then topint=ny-1
            if topint-bottomint+1 lt 5 then goto,skipthisnp
            xarr=findgen(topint-bottomint+1)
            
            params=[1.0,float(max(xarr)/2.0),npwidth/(2.0*SQRT(2.0*ALOG(2.0)))]
            ptemp=gaussian(xarr,params)
            ptemp=ptemp/total(ptemp)
            ;            print,ptemp
            bmep_mosdef_blind_extract,$
              npyexpect,$ ; yposition
              npwidth,$ ; width
              ny,ptemp,sciimg,var_img,f,ferr,fopt,fopterr,p
              
            bmep_mosdef_blind_save,savepath,maskname,filtername,slitname,$
              objnum,extrainfo1,extrainfo2,npyexpect,npwidth,0,$
              f,ferr,fopt,fopterr,p,min_width
              
            suffix='.'+ssi(objnum)
            if file_test(savepath+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits')  then begin
              data=readfits(savepath+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits',shdr,exten_no=5,/silent)
              cgplot,data,title=maskname+'.'+filtername+'.'+slitname+suffix+' y profile',ytitle='P'
              cgplot,(p/max(p))*max(data),color='red',/overplot
            endif
            
            
          endif;ct ne 0
          skipthisnp:
        endfor
      endif ;else print,'no np objects' ;ct gt 0
      
    ;      stop
      
    endif ; yexpect -1
  ;    endif ; NOREPEAT and file test
  ;    stop
  endfor ; n_ele objects
  
  
  theend:
  if norepeat eq 0 then ps_end
  cd,original_dir
  !p.multi=[0,1,1]
  print,'blind extraction took: ',round(systime(/seconds)-starttime),' seconds
  print,'end of best mosfire extraction program'
end





;TODO
;1.+ add in something special about the slit star
;  +enforce the width of the slit star to be the minimum width allowed
;  +calc the difference between expected and actual y position,
;  +use y pos difference to draw better line on 2d spectra
;  +case of 2 or more slit stars? avg the data? what do.
;  +external extra file to store this info?
;  +note:priority is -1 or 1 for stars.
;  +use this to find which are stars.
;  +fix naming for 2 objects in slit
;2.+ extract object by object, not filter by filter.
;  +the way i've done it is all of the objects in one filter
;  +then do the next filter.  This is not ideal
;  +should be obj by obj
;  +I need to think about how to implement this
;3. write program to re-extract objects.
;  +from wavel solution, calculate pixels to bin from
;  +bin the pixels, then fit a gaussian
;  +use expected position as initial guess
;  +check fwhm etc...
;  +comparison document/plots? (width, ypos, extraction...)
;  +what to do about stuff that was adjusted by hand?
;  +maybe have a flag on those?
;4. write blind extraction program
;  +get offset and fwhm from slit star.
;  +decide on a width? 2x? 1x?
;  +use the same profile for optimal weighting as slitstar
;  +do the extraction and save!
;  +if there is a secondary object in another slit,
;   extract at that expected position
;  Xmanually add in objects for blind extraction?
;5.+ for rotation or broad profiles, get all the flux.
;  +maybe edit fitting program to be a little bit wider.
;6.+write if width adjusted manually? or manual center'd?
;  +add in flag to not auto do this object?
;7.+ add wavel sol and all info to header. easy fix, i assume
;8.+make only essential save files. don't be lazy!!!
;9.+experiment with non-parametric width/center estimate
;  +just do this at the same time as normal gaussian fit
;  +only consider small region??
;10.+ to fit to a gaussian or not? that is the question?
;  +brian suggests fit to a  gaussian by default
;  +chisqr thresh to determine if gaussian shouldn't work?
;  +chi by eye?
;11. extract a sky spectrum?
;  -need to wait for mariska to add sky extension to 2d
;  -just extract at same position and width
;extras:
;  -plot of all of the data?
;  -compare box to optimal
;  -compare drp to mariska 2d/1d
;
;
;
;changes:
;added CRPIX1, CRTYPE1 to fits header
;added GAUSCHI to fits header. has chi sqr of gauss fit.
;fixed case of backwards l_bins and r_bins. now are sequential always.
;removed option to chose what filter to open
;instead, when opening an object all of the filters open
;title of window includes the filter
;
;added option to do "stars" as an object,
;this searches the SlitList.txt for something with priority 1 or -1
;throws fatal error if no stars found
;Added ISSTAR to fits header, 1 if star, 0 if not, -1 if something,
;bad happened
;00_starinfo.txt is created in the 2D folder if it does not exist at startup.
;00_starinfo.txt -maskname, filtername, objname, yexpect, yactual, width
;this file is read in for each object
;added isstar keyword to fits header. also used to determine if,
;main program should write to 00_starinfo.txt.
;main program now writes to 00_starinfo.txt if object is a star
;overwrites old info if star was already in list.
;
;wrapper averages star offset and width for all things with,
;matching maskname and filtername.
;
;from now on the yexpect keyword will be the one SHIFTED,
;from the star's center.  not the pure calculated one.
;
;reduced fitting width from 17 to 10 for gaussian.
;
;set default of pgauss to 1 instead of 0.
;
;30 ojects in 3 filters ~ 100 windows.   Too many! maybe only do,
;first 10 then second 10...
;^ fixed this by making user enter a min number and max number and,
;program does objects between those numbers.
;
;fixed .1 being added after primary of slit with multiple objects.
;
;added in "object number" or OBJNUM to fits header
;added objnum option in program by hitting 1,2,3,4,5.  This,
;sets the suffix and keyword.  This will allow for a secondary
;object to be extracted if there is no obvious primary.  This
;will also allow for if there are two very faint elines at different
;wavelengths one can bin separately for two different objects.
;
;got the extensions settled out...
;fixed the header to be on each extension.
;
;added comments describing each extension to fits hdr.
;
;added an auto extract flag that can be changed
;this is called AUTOEX in the fits header
;1 = yes, do the re-extraction automatically
;0 = no, do not do the re-extraction automatically.
;
;stars are skipped when doing all objects, they
;should be done already, so no sense in doing them again.
;
;confirmed for any object with a star, minimum width is applied
;
;fixed bug with min/max object num not actually working correctly.
;
;
;Issue with fits file extensions and headers NOT SOLVED completely,
;sent an email to mariska asking for help
;
;reduced width to fit gaussians to make chi sqr more reasonable (its now 8)
;
;switched the function gaussfit() for the function MPFITPEAK(),
;in the "p" calculating program.  gaussfit was causing underflow,
;errors for no apparent reason.
;
;fixed the fits headers w/ mariska's help
;
;changed y midpoint calculation back to integer math.
;
;removed duplicate "objnum" from fits header
;added comments to fits header
;
;removed the filter printing as it is no longer used. (no more select by filter)
;
;added width by hand and center by hand arrays in extraction code.
;added corresponding width by hand and center by hand keywords,
;in the fits headers (WBYHAND, CBYHAND)
;
;added calculation of the first two moments of profile for comparison with
;gaussian parameters.
;added MOM1 and MOM2 to fits header, which are these moments.
;
;fixed bug if there is a secondary object in the star frame, the yshift,
;and min width would be calculated from the secondary, not the primary.
;
;program now *automatically* detects if the monitor width is too small and,
;adjusts the click location accordingly.
;
;confirmed that two stars works as intended (yshift is avg, minw is smallest w)
;
;wbyhand keyword seems to work ok.
;
;fixed bug with two objects being extracted if one was not found w/ a gaussian fit
;fixed bug if objects were deleted some arrays should have also been deleted,
;but were not.
;
;created program called "bmep_mosdef_getinfo" which prints information about,
;EVERYTHING extracted. (except the blind stuffs)
;
;removed wavelength headers from the 5th extension.
;added spaces to comments in fits header
;
;updated blind extraction program to consider the shift from the stars
;updated blind extraction to do optimal extraction using star profile,
;5th extension is the weighting profile.
;
;program no longer outputs .txt files of the spectra.
;
;added isstar, yexpect, wbyhand, and objnum to blind header
;
;updated prerequisite comments above bmep program
;
;added in some default values for header keywords for blind extraction
;added 'BLIND' keyword to both headers (blind and normal)
;
;added ability to adjust the width of the window of the fit to automatically,
;determine the width of the profile (from a gaussian fit)
;
;added "norepeat" flag to beginning of blind extraction program.  This stops,
;the program from redoing objects that do not need to be blindly extracted.
;
;updated the crpix and other wavelength keywords for the blind extraction
;
;fixed bug with bad pixel mask not resetting for blind extraction
;fixed indexing error with optimal blind extraction
;added comparison plots to compare blind and normal extractions
;this automatically compares things blindly extracted vs non blind,
;IF THEY EXIST!
;
;alphabetized subroutines
;
;added FORWARD_FUNCTION statements to start of almost every subroutine due to,
;conflict created when I alphabetized everything.
;
;bkgnd subtracted image not plotted if state.savetext eq 0. this,
;stops creation of duplicate plot for mosfire data (since its already bkgnd subbed)
;
;swapped the sort() function in bmep_mosdef_getinfo to bsort() due to IDL lame-ness
;also fixed sorting issue in bmep_mosdef_getinfo
;
;fixed bug in blind extraction if star slit width was smaller than actual slit making the profile get messed up.
;
;put some sections of the blind extraction into their own subroutines.
;blind program creates 00_np_info.txt (non primary) that has information,
;about the non-primary objects for the purposes of blind extraction.
;
;made the 'x' key interactive instead of typing in the wavelengths.
;double tapping this key or having the first position to the right
;of the second position will cause the range to be reset
;
;The 'Z' key now guesses the redshift of an object by taking into consideration
;the z of the field, the filter, and the location of the brightest line extracted.
;This doesn't work if you extracted only wide wavelength ranges.  There must be a
;line to consider.
;
;The places where you clicked to bin to get the profile are shown in red
;on the extracted profile.  These are numbered (though its hard to see)
;
;The 'X' key works kind of like the 'g' key in iraf.  This fits gaussian to
;the data inbetween the two points that your mouse is hovering over. if
;its a good fit you have the option to select which line this corresponds to or
;input a wavelength of a line.  Then this information is saved to the headers of the
;saved fits files (still working on that last bit)
;
;changed programs with "mariska" in their name to have "mosdef" instead
;
;jan 2014
;
;added calculation of slitloss programs (bmep_calc_slitloss*)
;added more parallel arrays with information about the fit and slitloss
;parallel arrays are saved to header
;
;now the program ONLY SAVES .FITS files, including all info in hdr.
;program MUST HAVE maskname, slitname, filtername in extrainfo1,2,3
;
;now saves a 6th extension, the error in the yprofile...
;
;added sigmastar to 00_starinfo.txt. changed center value to actually calculated center.
;Added support for minw to be a float.  It used to be rounded at a bunch of steps.
;
;march 2014
;
;Revamped y range to be more reasonable if zoomed in. (threshold is 300 pixels)
;
;renamed var_img to std_img to  be more clear. The var_img is then calculated from this
;
;renamed starinfo.txt to be 00_starinfo.txt
;
;may 2014
;



pro bmep_mosdef,path_to_output=path_to_output
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  !except=2 ;see division by zero errors instantly.
  astrolib
  
  ;set output to what is in the envoirnment variable
  path_to_output=getenv('BMEP_MOSDEF_2D')
  if x ne '' then path_to_output=x
  ;default if no env found
  if ~keyword_set(path_to_output) then path_to_output='~/mosfire/output/idl_output/2D/' ; no trailing slash
  
  ;ensure that there is a '/' at the end of the path.
  if strmid(path_to_output,strlen(path_to_output)-1) ne '/' then path_to_output=path_to_output+'/'
  cd,path_to_output,current=original_dir
  
  ;clear away all windows!!
  close,/all
  ireset,/no_prompt
  
  ;get where to save output of extraction program
  x=getenv('BMEP_MOSDEF_1D')
  if x ne '' then savepath=x else begin
    savepath=path_to_output+maskname+'/1d_extracted/'
    ;create folder to extract to if it doesn't exist.
    if ~bmep_DIR_EXIST(savepath) then spawn,'mkdir 1d_extracted'
    endelse
  ;ensure that there is a '/' at the end of the path.
  if strmid(savepath,strlen(savepath)-1) ne '/' then savepath=savepath+'/'
  
  ;create a 00_starinfo.txt if not exist
  if ~file_test(savepath+'00_starinfo.txt') then begin
    forprint,textout=savepath+'00_starinfo.txt',['maskname '],$
      ['filtername '],['objname '],[99.99],[99.99],[99.99],[99.99],/nocomment
  endif
  
  
  ;parse folder into different masks!!
  spawn,'ls *.2d.fits > tempfiles.txt'
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt'
  
  ;parse names
  masks=[]
  filters=[]
  objects=[]
  for i=0,n_elements(filenames)-1 do begin
    substrings=strsplit(filenames[i],'.',/extract)
    if n_elements(substrings) eq 5 then begin
      masks=[masks,substrings[0]]
      filters=[filters,substrings[1]]
      objects=[objects,substrings[2]]
    endif ;else message,'ERROR something without a proper name is in the folder.'
  end
  
  ; user  choose mask
  masks_no_duplicates=masks[rem_dup(masks)]
  forprint,indgen(n_elements(masks_no_duplicates)),' '+masks_no_duplicates
  print,'which number'
  read,choice
  if choice lt 0 then goto,theend
  if choice ge n_elements(masks_no_duplicates) then goto,theend
  maskname=masks_no_duplicates[choice]
  print,'mask: ',maskname
  
  ;remove everything except that mask from "database"
  index=where(masks eq maskname,ct)
  if ct eq 0 then message,'error' ; check for strange errors
  masks=masks[index]
  filters=filters[index]
  objects=objects[index]
  
  
;  ; user  choose filter.
;  filters_no_duplicates=filters[rem_dup(filters)]
;  if n_elements(filters_no_duplicates) ne 1 then begin
;    ;  forprint,indgen(n_elements(filters_no_duplicates)),' '+filters_no_duplicates
;    ;  print,'which number'
;    ;  read,choice
;    choice=0
;    if choice lt 0 then goto,theend
;    if choice ge n_elements(filenames) then goto,theend
;  endif else choice=0
;  filtername=filters_no_duplicates[choice]
;  
;  ;remove everything except that filter from "database"
;  index=where(filters eq filtername,ct)
;  if ct eq 0 then message,'error' ; check for strange errors
;  masks=masks[index]
;  filters=filters[index]
;  objects=objects[index]
  
  
  ;user choose object.
  forprint,indgen(n_elements(objects)),' '+objects
  print,n_elements(objects),' all files'
  print,n_elements(objects)+1,' all stars'
  print,'which number'
  choice=0
  read,choice
  if choice lt 0 then goto,theend
  if choice ge n_elements(objects)+2 then goto,theend
  ;do all objects
  if choice eq n_elements(objects) then begin
    choicearr=indgen(n_elements(objects))
    minobjnum=0
    maxobjnum=n_elements(objects)-1
    print,'min object num?'
    read,minobjnum
    print,'max object num?'
    read,maxobjnum
    index=where(choicearr ge minobjnum and choicearr le maxobjnum,/null)
    choicearr=choicearr[index]
  endif else choicearr=[choice]
  
  ;FIND ALL THE STARS!!!
  if choice eq n_elements(objects)+1 then begin
    slitlistfile=path_to_output+'00mask_info/'+maskname+'_SlitList.txt'
    if ~file_test(slitlistfile) then message,'SLITLIST FILE MISSING... '+SLITLISTFILE
    readcol,slitlistfile,slitnamearr,priorityarr,offsetarr,format='X,X,X,X,X,X,X,X,X,A,F,F,X,X,X,X,X,X',/silent
    star_index=where(abs(priorityarr) eq 1.0,ct)
    if ct eq 0 then message,'no stars???? check '+slitlistfile+' for an object of priority 1 or -1'
    print,slitnamearr[star_index],' are stars'
    choicearr=[]
    for k=0,n_elements(star_index) - 1 do begin
      index=where(sss(objects) eq sss(slitnamearr[star_index[k]]))
      choicearr=[choicearr,index[0]]
    endfor
  endif ; chose to do all stars.
  
  for i=0,n_elements(choicearr)-1 do begin
    for j=0,n_elements(filters_no_duplicates)-1 do begin
      slitname=objects[choicearr[i]]
      filtername=filters_no_duplicates[j]
      filename=maskname+'.'+filtername+'.'+slitname+'.2d.fits'
      print,filename
      
      if ~file_test(filename) then begin
        print,'WARNING: NO IMAGE FOUND, ',FILENAME
        print,'WARNING: NO IMAGE FOUND, ',FILENAME
        print,'WARNING: NO IMAGE FOUND, ',FILENAME
        goto, no_do_image
      ENDIF
      
      ;read in files
      sciimg=readfits(filename,shdr, /SILENT,exten_no=1)
      sciimg=double(sciimg)
      
      index=where(finite(sciimg) eq 0,/null)
      sciimg[index]=0.0
      
      ;calculate variance image
      noise_img=readfits(filename, /SILENT,exten_no=3)
      ;clean image
      index=where(noise_img eq 99,/null)
      sciimg[index]=0.0
      noise_img[index]=0.0
      index=where(finite(noise_img) eq 0,/null)
      sciimg[index]=0.0
      noise_img[index]=0.0
      noise_img=double(noise_img)
      var_img=noise_img*noise_img
      
      ;calculate wavel.
      wavel=(sxpar(shdr,'CRVAL1')+findgen(n_elements(sciimg[*,0]))*sxpar(shdr,'CDELT1'))
      
      ;combine images and make snr image
      ny=n_elements(sciimg[0,*])
      nx=n_elements(sciimg[*,0])
      
      ;snr image
      index=where(var_img eq 0,ct)
      if ct ne 0 then var_img[index]=9801.0
      snrimg=abs(sciimg/sqrt(var_img))
      snr2sigCut=snrimg
      index=where(snrimg lt 2.,ct)
      if ct gt 0 then snr2sigCut[index]=0.0

      
      del=10.
      botpercent=del
      toppercent=100.-del
      
      kernel = GAUSSIAN_FUNCTION([500,5], WIDTH=15, MAXIMUM=255,/double)
      kernel=kernel-avg(kernel)
      kernel=[[0,-1,0],$
              [-1,5,-1],$
              [0,-1,0]]      
      kernel=[[-1,-1,-1],$
              [-1,8,-1],$
              [-1,-1,-1]]     
      kernel=findgen(11,11)
      for ii=0,n_elements(kernel[*,0])-1 do $
        for jj=0,n_elements(kernel[*,0])-1 do $
        kernel[ii,jj]=(-1.0)*bmep_LoG(ii-n_elements(kernel[*,0])/2,jj-n_elements(kernel[*,0])/2,3) 
;      kernel=[[ 0,-1,-1,-1, 0],$
;              [-1, 3, 3, 3,-1],$
;              [-1, 3, 5, 3,-1],$
;              [-1, 3, 3, 3,-1],$
;              [ 0,-1,-1,-1, 0]]      
;      kernel=[[ 0,-1,-1,-1, 0],$
;              [-1, 3, 5, 3,-1],$
;              [-1, 5, 7, 5,-1],$
;              [-1, 3, 5, 3,-1],$
;              [ 0,-1,-1,-1, 0]]
      kernel=kernel-avg(kernel)
      
;      print,kernel
      conv_img=convol((sciimg),kernel,/edge_zero,/normalize)
      for ii=0,n_elements(conv_img[*,0])-1 do conv_img[ii,*]=conv_img[ii,*]-aVG(conv_img[ii,*])

      big_img=findgen(nx,(ny*4))
      big_img[*,ny*0:ny*1-1]=bytscl(sciimg,top=255,/NAN,$
        min=bmep_percent_cut(sciimg,botpercent),max=bmep_percent_cut(sciimg,toppercent))   ;science img
      big_img[*,ny*1:ny*2-1]=bytscl(conv_img,top=255,/NAN,$
        min=bmep_percent_cut(conv_img,00),max=bmep_percent_cut(conv_img,100))  ;var img 0,15
      big_img[*,ny*2:ny*3-1]=bytscl(snrimg,top=255,/NAN,$;min=98.5,max=99.0)
        min=bmep_percent_cut(snrimg,botpercent),max=bmep_percent_cut(snrimg,90.0))   ;snr img
      big_img[*,ny*3:ny*4-1]=bytscl(snr2sigCut,top=255,/NAN,$
        min=0.0,max=3.5)   ;snr img
        
        
      ;calculate where object SHOULD be!
      slitlistfile=path_to_output+'00mask_info/'+maskname+'_SlitList.txt'
      if ~file_test(slitlistfile) then message,'SLITLIST FILE MISSING... '+SLITLISTFILE
      
      ;1 2 17  45.61 -5  10  9.37  0.70  7.01  23763 800.00  -0.38 2 17  45.59 -5  10  9.65
      readcol,slitlistfile,slitnamearr,priorityarr,offsetarr,format='X,X,X,X,X,X,X,X,X,A,F,F,X,X,X,X,X,X',/silent
      index=where(sss(slitnamearr) eq sss(slitname),ct)
      yexpect=-1
      isstar=-1
      minwidth=-1
      yshift=0.0
      if ct eq 1 then begin
        pixscale=0.1799
        midpoint=ny/2
        yexpect=midpoint+offsetarr[index[0]]/pixscale
        ;check if object is a star
        if abs(priorityarr[index]) eq 1 then isstar=1 else begin
          isstar=0
          readcol,savepath+'00_starinfo.txt',maskstar,$
            filtstar,objstar,yexpect_star,yactual_star,widthstar,sigmastar,/silent,format='A,A,A,F,F,F,F'
          index=where(maskstar eq maskname and filtstar eq filtername,ct)
          if ct ne 0 then begin
            print,ct,' number of stars found for ',maskname,' ',filtername
            yshift=avg(yexpect_star[index] - yactual_star[index])
            minwidth=min(widthstar[index])
          endif else print,'no matching stars found ',maskname,' ',filtername
        endelse ;if isstar
        
        if isstar eq 1 and choice eq n_elements(objects) then goto, no_do_image
        print,'slitname, yexpect, midpoint, yshift, minwidth'
        print,slitname, yexpect, midpoint, yshift, minwidth
        
        yexpect=yexpect-yshift
        PRINT,'new yexpect:',yexpect
      endif else print,'no object found in the slitlist file?!?!?!?'
      
      ;draw white line
      big_img[*,yexpect]=255
      
      highval=max(big_img)
      lowval=min(big_img)
      
      
      ;calculate info to add to hdr (include zguess placeholders)
      extrainfo1=[$
        'CRVAL1',$
        'CDELT1',$
        'CRPIX1',$
        'CTYPE1',$
        'EXPTIME',$
        'FILNM',$
        'MSKNM',$
        'FILTNM',$
        'SLITNM',$
        'ISSTAR',$
        'MINW',$
        'YEXPECT'$
        ]
        
      extrainfo2=[$
        string(sxpar(shdr,'CRVAL1')),$
        string(sxpar(shdr,'CDELT1')),$
        string(sxpar(shdr,'CRPIX1')),$
        'LINEAR',$
        string(sxpar(shdr,'EXPTIME')),$
        filename,$
        maskname,$
        filtername,$
        slitname, $
        ssi(isstar), $
        ssf(minwidth), $
        ssf(yexpect) $
        ]
        
      ;comments
      extrainfo3=[$
        ' ',$
        ' ',$
        ' ',$
        ' ',$
        ' total exposure time (seconds)',$
        ' name of file',$
        ' name of mask',$
        ' name of filter',$
        ' name of slit', $
        ' Flag if is a star', $
        ' minimum width (-1 default)', $
        ' expected y position' $
        ]
        
      bmep_display_image,big_img,sciimg,var_img,highval,lowval,slitname,filtername,wavel,savepath,$
        revisevar=0,extrainfo1=extrainfo1,extrainfo2=extrainfo2,extrainfo3=extrainfo3,savetext=0
      
;      cghistoplot,snrimg,xr=[0,15],binsize=0.2,ytitle='SIGNAL TO NOISE RATIO',TITLE=MASKNAME+FILTERNAME+SLITNAME
;      stop
      
      no_do_image:
    endfor ; filters_no_duplicates
  endfor ; choicearr
  
  theend:
  
  cd,original_dir
  print,'end of best mosfire extraction program'
end

FUNCTION bmep_MouseDown, oWin, x, y, iButton, KeyMods, nClicks
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  state = oWin.UVALUE
  state.x0 = x
  state.y0 = y
  print,'mouse down at',x,y
  state.buttonDown = 1
  oWin.UVALUE=state
  RETURN, 0 ; Skip default event handling
END


FUNCTION bmep_MouseUp, oWin, x, y, iButton
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  dimensions = GET_SCREEN_SIZE()
  monitor_width=dimensions[0]-10
  state = oWin.uvalue
  IF (~state.buttonDown) THEN RETURN, 0
  x0 = state.x0
  y0 = state.y0
  print,'mouse up at',x,y
  print,x0,' to ',x
  
  if n_elements(state.data[*,0]) gt monitor_width then begin
    print,'monitor not wide enough ', (float(n_elements(state.data[*,0]))/float(monitor_width))
    print,'trying to account for this...'
    x0=x0* (float(n_elements(state.data[*,0]))/float(monitor_width))
    x =x * (float(n_elements(state.data[*,0]))/float(monitor_width))
    print,x0,' to ',x
  endif
  
  if x0 eq x then x=x+1
  
  ;check boundries
  if x  gt n_elements(state.data[*,0])-1 then  x=n_elements(state.data[*,0])-1
  if x0 gt n_elements(state.data[*,0])-1 then x0=n_elements(state.data[*,0])-1
  if x  lt 0 then  x=0
  if x0 lt 0 then x0=0
  
  
  
  ;plot a cross section of area between mouse clicks
  for i=0,n_elements(state.ydata)-1 do begin
    ;add each individual point
    if x0 ne x then begin
      if x0 lt x then index=indgen(abs(x-x0))+x0 else index=indgen(abs(x-x0))+x
      
      if state.cont_mode eq 1 then begin
        noiseguess=[]
        for j=0,n_elements(index)-1 do noiseguess=[noiseguess,total(state.var_img[index[j],*])]
        goodpix=where(noiseguess lt median(noiseguess))
        badpix=where(noiseguess ge median(noiseguess))
        
        ;          if i eq 0 then begin
        ;            print,'including columns:'
        ;            index2=index
        ;            index2[badpix]=0
        ;            forprint,index2[goodpix],noiseguess[goodpix]
        ;            print,n_elements(goodpix),n_elements(badpix)
        ;            endif
        
        index=index[goodpix]
      endif
      
      
      for j=0,n_elements(index)-1 do state.ydata[i]=state.ydata[i]+state.data[index[j],i]
      for j=0,n_elements(index)-1 do state.ydataerr[i]=state.ydataerr[i]+state.var_img[index[j],i]
    endif
  endfor
  
  ;;normalize
  ;state.ydata=state.ydata/max(state.ydata)
  
  plot,state.ydata,/ynozero
  oplot,replicate(0.0,n_elements(state.data)),color=255.+255.*255.
  index=where(state.extrainfo1 eq 'YEXPECT',ct)
  if ct eq 1 then oplot,[float(state.extrainfo2[index]),float(state.extrainfo2[index])],minmax(state.ydata),color=255.+255.*255.
  
  ;add in the collapsed area to header info
  if state.n_bins ge 10 then print,'too many bins' else begin
    if x0 le x then begin
      state.l_bins[state.n_bins]=x0
      state.r_bins[state.n_bins]=x
    endif else begin
      state.l_bins[state.n_bins]=x0
      state.r_bins[state.n_bins]=x
    endelse
    state.n_bins=state.n_bins+1
  endelse
  
  
;stats_mode stuff
;ny=n_elements(state.data[0,*])
;  print,'stats from ',x0,x
;  print,'stats from ',y0,y
;if y gt y0 and y0 lt ny then begin
;  print,'stats from ',x0,x
;  print,'stats from ',y0,y
;  small_img=   state.data[x0:x,y0:y]
;  small_var=state.var_img[x0:x,y0:y]
;  print,avg(small_img),' avg signal'
;  print,avg(small_var),' avg var'
;  print,avg(sqrt(small_var)),' avg std'
;  endif
  
  
  
  
  
  oWin.uvalue=state
  
  ; Clear the current selections
  oSelect = oWin.GetSelect()
  FOREACH oVis, oSelect do oVis.Select, /UNSELECT
  RETURN, 0 ; Skip default event handling
  
END


;return the value of an array at the nth percentile.
function bmep_percent_cut,array,percent,above=above,below=below
  ;check inputs
  if percent le 0 then return,min(array)
  if percent ge 100 then return,max(array)
  
  ;temp var to not mess with input
  x=array
  
  ;eleminate bad pixels
  index=where(finite(x) ne 0,ct)
  if ct eq 0 then return,-999999
  x=x[index]
  
  ;sort data and give a value at the nth percentile
  percentfrac=double(percent)/100.0
  nele=double(n_elements(x))
  index=sort(x)
  sarr=x[index];ascending order
  return,sarr[long(nele*percentfrac)]
end



;plot the statistics of an area selected by a user.
pro bmep_plot_stats,state,leftstatspos,rightstatspos,f,fopt,statshift,ferr,fopterr

  ;find area of interest.
  index=where(state.wavel ge leftstatspos and state.wavel le rightstatspos,ct)
  
  ;make sure there are enough stats points.
  if ct le 3 then print,'no flux between cursor marks' else begin
    momentresult=moment(f[index],sdev=std_dev)
    momentresultopt=moment(fopt[index],sdev=std_devopt)
    statspacing=0.05
    statbeginning=0.90
    
    ;put stats on the plot!
    xyouts,0.05,statbeginning-statspacing*(-1.0)-statshift,$
      '           boxcar        optimal',/normal,charsize=2
    xyouts,0.05,statbeginning-statspacing*0.0-statshift,$
      'mean    :'+ssf(momentresult[0])+'... '+ssf(momentresultopt[0]),/normal,charsize=2
    xyouts,0.05,statbeginning-statspacing*1.0-statshift,$
      'std dev :'+ssf(std_dev)+'... '+ssf(std_devopt),/normal ,charsize=2
    xyouts,0.05,statbeginning-statspacing*2.0-statshift,$
      'SNR sigma:'+ssf(momentresult[0]/std_dev)+'... '+$
      ssf(momentresultopt[0]/std_devopt),/normal,charsize=2
    xyouts,0.05,statbeginning-statspacing*3.0-statshift,$
      'SNR ebar:'+ssf(momentresult[0]/avg(ferr[index]))+'... '+$
      ssf(momentresultopt[0]/avg(fopterr[index])),/normal,charsize=2
    statshift=statshift-statspacing*3.0
    
    ;plot vertical lines
    oplot_vert,leftstatspos,minmax(f),3
    oplot_vert,rightstatspos,minmax(f),3
  endelse ; ct gt 3
end








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






function bmep_sigma_clip,arr,wind,N_sigma=N_sigma
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  n_replaced=0
  if ~keyword_set(N_sigma) then N_sigma=2
  for i=0,n_elements(arr)-wind,wind-wind/2 do begin
    range=indgen(wind)+i
    sdev=stddev(arr[range])
    smallarr=arr[range]-median(arr[range])
    index=where(abs(smallarr) gt sdev*N_sigma,ct)
    if ct gt 0 and sdev ne 0 then arr[range[index]]=median(arr[range])
    n_replaced+=ct
  endfor
  print,n_replaced
  return,arr
end


pro bmep_test
print,0
;test idl version number.
IF (Float(!Version.Release) lt 8.1) THEN message,'idl ver must be 8.1 or greater'
print,1
x=1.11
y=ssf(x);test small programs 
y=ssi(x);test small programs 
y=sss(x);test small programs 
y=nearest_to(findgen(100),49.5);test small programs 
print,2
astrolib;test astrolib
x=minmax(findgen(100));test astrolib
print,3
plot,findgen(100);check if plotting works
oplot_vert,20,[0,100]
print,4

x=image(findgen(100,100))
x.close
print,5
print,'BMEP_MOSDEF_2D is ',getenv('BMEP_MOSDEF_2D')
print,'BMEP_MOSDEF_1D is ',getenv('BMEP_MOSDEF_1D')
print,'BMEP_MOSDEF_MASKS is ',getenv('BMEP_MOSDEF_MASKS')
print,'BMEP_MOSFIRE_DRP_2D is ',getenv('BMEP_MOSFIRE_DRP_2D')
print,'BMEP_MOSFIRE_DRP_1D is ',getenv('BMEP_MOSFIRE_DRP_1D')
print,'BMEP_MOSFIRE_DRP_MASKS is ',getenv('BMEP_MOSFIRE_DRP_MASKS')
print,'BMEP_DROPBOX_PATH is ',getenv('BMEP_DROPBOX_PATH')
print,6

print,7

print,8

print,9

print,10



print,'test over'
end


;view 1d mosfire files.
pro bmep_view_1d,path_to_dropbox=path_to_dropbox,path_to_output=path_to_output
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
    
  ;defaults
  astrolib
  if not keyword_set(path_to_dropbox) then path_to_dropbox='~/Dropbox/'
  if not keyword_set(path_to_output) then path_to_output='~/mosfire/output/'
  cd,path_to_output,current=original_dir
  
  ;clear away all windows!!
  close,/all
  ireset,/no_prompt
  
  ;find folders of masks..
  spawn,'pwd'
  spawn,'ls > tempfiles.txt'
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt
  
  ;print out folders and ask user which to use
  forprint,indgen(n_elements(filenames)),replicate(' ',n_elements(filenames)),filenames
  print,'which number'
  read,choice
  
  ;if choice was bad then bail. Yeah. goto statements. DEAL WITH IT.
  if choice lt 0 then goto,theend
  if choice ge n_elements(filenames) then goto,theend
  
  ;foldername is same as the mask name
  maskname=filenames[choice]
  print,'mask: ',maskname
  cd,maskname
  savepath=path_to_output+maskname+'/1d_extracted/'
  if not bmep_DIR_EXIST(savepath) then message,'no 1d extracted for this mask'
  cd,savepath
  
  ;read in extracted
  spawn,'ls *.sav > tempfiles.txt'
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt
  
  ;get which file to view
  forprint,indgen(n_elements(filenames)),replicate(' ',n_elements(filenames)),filenames
  print,'which number'
  read,choice
  
  ;if choice was bad then bail. Yeah. goto statements. DEAL WITH IT.
  if choice lt 0 then goto,theend
  if choice ge n_elements(filenames) then goto,theend
  
  restore,filenames[choice]
  
  ;view 2d image
  newimg=bytscl(img2d,top=255,/NAN,$
    min=bmep_percent_cut(img2d,10),max=bmep_percent_cut(img2d,90))
  im=image(newimg,dimensions=[n_elements(newimg[*,0]),n_elements(newimg[0,*])],$
    window_title=slitname,margin=[0,0,0,0],max_value=255,min_value=0)
    
  ;view 1d using iraf_plot.
  iraf_plot,wavel,FLUX_OPT,returnval,title=SLITNAME,$
    plot_title=SLITNAME, path_to_dropbox=path_to_dropbox
    
  ;close window and clean up
  im.window.close
  theend:
  cd,original_dir
  print,'end of bmep_view_1d'
end


;rereduce all data for lris.
pro bmep_wrap_rereduce_optimal,path_to_output=path_to_output
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  astrolib
  
  if not keyword_set(path_to_output) then path_to_output='~/LRIS/'
  cd,path_to_output,current=original_dir
  
  ;find folders of masks..
  spawn,'ls > tempfiles.txt'
  readcol,'tempfiles.txt',foldernames,format='A',/silent
  spawn,'rm tempfiles.txt
  
  print,'foldernames'
  print,foldernames
  print
  
  for i=0,n_elements(foldernames)-1 do begin
    extracted_folder=foldernames[i]+'/1d_extracted'
    
    ;if the folder is real, extract!
    if bmep_dir_exist(extracted_folder) then begin
      cd,extracted_folder
      filenames=[]
      spawn,'ls *.sav > tempfiles.txt'
      readcol,'tempfiles.txt',filenames,format='A',/silent,count=count
      spawn,'rm tempfiles.txt
      
      ;if any .save files are found, rereduce them!
      for j=0,n_elements(filenames)-1 do begin
        print
        print,'bmep_rereduce_optimal - ',extracted_folder+filenames[j]
        bmep_rereduce_optimal,filenames[j]
      endfor;num of files
    endif ; dir exists
    cd,path_to_output
  endfor;foldernames
  print,'end of bmep rereduce wrapper'
end




;============================================================================
;
;
;Best Mosfire Extraction Program!! (also made to work with LRIS data.)
;billfreeman44@yahoo.com
;
;============================================================================
;this program is supposed to make the reduction of spectroscopic data easier.
;needs idl 8.1 or greater!!!!!!!!!!!!!!11
;needs to read the output folders for MOSFIRE
;note: it deletes all windows from IDL at startup. (plots and images)
;compile twice before running.
;needs:
; -astronomy library
; -mpfit library
; -coyote library
; -programs i've written:
;    oplot_vert.pro
;    ssf.pro
;    ssi.pro
;    sss.pro
;    nearest_to.pro
;for bmep mosdef:
;  2d spectra should be in '~/mosfire/output/idl_output/2D'
;
;  mask in '~/mosfire/output/idl_output/2D/mask_info'
;  the *slitlist.txt files should have no S in them,
;  by the stars in the slitname columns.
;
;  there should be '~/mosfire/output/idl_output/2D/mask_info/2012B'
;  and '~/mosfire/output/idl_output/2D/mask_info/2013A' folders.
;
;
;UPDATES:
; added a "gzending" keyword because apparently the DRP is not saving as
; a .fits.gz file its saving as a .fits without the .gz...  if your files have a
; .gz at the end then set this keyword.
; 
; added a 'ivarending' keyword because the DRP seemed to stop 
; using the IVAR and is now using std or something
;
;test
;
pro bmep,path_to_dropbox=path_to_dropbox,path_to_output=path_to_output,gzending=gzending,ivarending=ivarending
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  !except=2

  
  astrolib
  if ~keyword_set(gzending) then gzending=0
  if ~keyword_set(ivarending) then ivarending=0
  if ~keyword_set(path_to_dropbox) then path_to_dropbox='~/Dropbox/'
  
  ;set output to what is in the envoirnment variable
  x=getenv('BMEP_MOSFIRE_DRP_2D')
  if x ne '' then path_to_output=x
  ;default if no env found
  if ~keyword_set(path_to_output) then path_to_output='~/mosfire/output/'
  
  ;ensure that there is a '/' at the end of the path.
  if strmid(path_to_output,strlen(path_to_output)-1) ne '/' then path_to_output=path_to_output+'/'
  print,'the output 2D path is',path_to_output
  
  cd,path_to_output,current=original_dir
  
  ;clear away all windows!!
  close,/all
  ireset,/no_prompt
  
  
  ;find folders of masks..
  spawn,'pwd'
  spawn,'ls > tempfiles.txt'
  readcol,'tempfiles.txt',filenames,format='A',/silent
  spawn,'rm tempfiles.txt
  forprint,indgen(n_elements(filenames)),replicate(' ',n_elements(filenames)),filenames
  print,'which number'
  read,choice
  if choice lt 0 then goto,theend
  if choice ge n_elements(filenames) then goto,theend
  
  
  ;assume foldername is same as the mask name
  maskname=filenames[choice]
  print,'mask: ',maskname
  cd,maskname
  
  x=getenv('BMEP_MOSFIRE_DRP_1D')
  if x ne '' then savepath=x else begin
    savepath=path_to_output+maskname+'/1d_extracted/'
    ;create folder to extract to if it doesn't exist.
    if ~bmep_DIR_EXIST(savepath) then spawn,'mkdir 1d_extracted'
    endelse
    
  ;create a 00_starinfo.txt if not exist
  if ~file_test(path_to_output+'00_starinfo.txt') then begin
    forprint,textout=path_to_output+'00_starinfo.txt',['maskname '],$
      ['filtername '],['objname '],[99.99],[99.99],[99.99],[99.99],/nocomment
  endif  
  
  ;check if fits files exist and create filenames array
  if gzending eq 0 then spawn,'ls '+maskname+'_?_*eps.fits > xx__spec_list.txt' $
  else spawn,'ls '+maskname+'_?_*eps.fits.gz > xx__spec_list.txt'
  readcol,'xx__spec_list.txt',filenames,format='A',/silent
  spawn,'rm xx__spec_list.txt
  if n_elements(filenames) eq 0 then message,'no fits files found... probably bad fliter choice'
  
  ;get the names of the slits
  slitnames=[]
  for i=0,n_elements(filenames)-1 do $
    slitnames=[slitnames,bmep_get_slitname(filenames[i],maskname,/eps,gzending=gzending)]
  slitnames_nodup=slitnames[rem_dup(slitnames)]

  
  
  ;get which slit to do.
  forprint,indgen(n_elements(slitnames_nodup)),replicate('. ',n_elements(slitnames_nodup)),slitnames_nodup
  print,n_elements(filenames),' all files'
  print,'which slit?'
  print,'The blank slit is the full 2d image... dont choose this one'
  print
  choice=0
  read,choice
  if choice lt 0 then goto,theend
  if choice ge n_elements(filenames)+1 then goto,theend
  
  if choice eq n_elements(filenames) then choicearr=indgen(n_elements(filenames)) else choicearr=where(slitnames eq slitnames_nodup[choice])
  
  for i=0,n_elements(choicearr)-1 do begin
    slitname=slitnames[choicearr[i]]
    epsfile=filenames[choicearr[i]]
    ;check if the files exist
    if gzending eq 0 then $
      if ivarending eq 0 then $
        stdfile=strmid(epsfile,strlen(epsfile)-strlen('_eps.fits'))+'_std.fits' $
        else $
          stdfile=strmid(epsfile,strlen(epsfile)-strlen('_eps.fits'))+'_ivar.fits' $
      else $
      stdfile=strmid(epsfile,strlen(epsfile)-strlen('_eps.fits.gz'))+'_ivar.fits.gz'
    ;check if these files actually exist
    if ~file_test(epsfile) then print,'WARNING!!!! file '+epsfile+' does not exist'
    if ~file_test(stdfile) then print,'WARNING!!!! file'+stdfile+' does not exist'
    
    if file_test(epsfile) and file_test(stdfile) and slitname ne '' then begin
      ;read in files
      sciimg=readfits(epsfile,shdr, /SILENT)
      sciimg=double(sciimg)
      std_img=readfits(stdfile, /SILENT)
      std_img=double(abs(std_img))
      
      ny=n_elements(sciimg[0,*])
      nx=n_elements(sciimg[*,0])
      
      ;fix the image if its an ivar image
      if ivarending eq 1 then begin 
        index=where(std_img ne 0.0,ct)
        std_img[index]=sqrt(1.0/std_img[index])
        endif  
      
      ;calculate variance image
      index=where(std_img ne 0.0,ct)
      var_img=replicate(0.0,nx,ny)
      var_img[index]=std_img[index]*std_img[index]
      
      ;clean variance image of nan or inf values
      ind=where(finite(var_img) eq 0,ct)
      if ct ne 0 then var_img[ind]=0.0
      
      
      ;FOR LONGSLIT CORRECTION TO MAKE THEM NOT HUGE IMAGES
;      sciimg=sciimg[*,977:1095]
;      std_img=std_img[*,977:1095]
;      var_img=var_img[*,977:1095]
;      ny=n_elements(sciimg[0,*])
;      nx=n_elements(sciimg[*,0])
      
      ;calculate wavel. CTYPE1, CRVAL1, and CDELT1
      refpix = sxpar(shdr,'CRPIX1')
      lam0   = sxpar(shdr,'CRVAL1')
;      print,'starting wavelength ', lam0
      delta  = sxpar(shdr,'CDELT1')
      if delta lt 1.1 then delta  = sxpar(shdr,'CD1_1')
      wavel=lam0 + (findgen(n_elements(sciimg[*,0]))+refpix-1.0) * delta
      
;      print,'CRPIX1: ',refpix
;      print,'CRVAL1: ',lam0
;      print,'CDELT1: ',sxpar(shdr,'CDELT1')
;      print,'CD1_1:  ' ,sxpar(shdr,'CD1_1')
;      PRINT,'DELTA:  ',delta
      
      ;snr image
      index=where(std_img ne 0)
      snrimg=sciimg/std_img[index]
      snr2sigCut=snrimg
      index=where(snr2sigCut lt 2.0,/null)
      snr2sigCut[index]=0.0
      
      ;mess with these to change how an image is viewed. (scaling)
      botpercent=10.0
      toppercent=90.0
      
      big_img=findgen(nx,(ny*4))
      big_img[*,ny*0:ny*1-1]=bytscl(sciimg,top=255,/NAN,$
        min=bmep_percent_cut(sciimg,botpercent),max=bmep_percent_cut(sciimg,toppercent))   ;science img
      big_img[*,ny*1:ny*2-1]=bytscl(ivarimg,top=255,/NAN,$
        min=bmep_percent_cut(ivarimg,botpercent),max=bmep_percent_cut(ivarimg,toppercent))  ;ivar img
      big_img[*,ny*2:ny*3-1]=bytscl(sciimg,top=255,/NAN,$
        min=bmep_percent_cut(abs(sciimg),5.0),max=bmep_percent_cut(abs(sciimg),95.0))   ;snr img
      big_img[*,ny*3:ny*4-1]=bytscl(snr2sigCut,top=255,/NAN,$
        min=2.0,max=3.0)   ;snr img
        
      ;resize big_img if really big in y direction...
      if ny*4 gt 1000 then begin
        big_img=bytscl(sciimg,top=255,/NAN,$
          min=bmep_percent_cut(sciimg,botpercent),max=bmep_percent_cut(sciimg,toppercent))
      endif
      
      ;find yexpect
      slitlistfile='/Users/bill/mosfire/output/idl_output/2D/mask_info/non_mosdef/masks/'+$
        maskname+'/'+maskname+'_SlitList.txt'
      yexpect=-1
      if file_test(slitlistfile) then begin
        readcol,slitlistfile,slitnamearr,priorityarr,offsetarr,format='X,X,X,X,X,X,X,X,X,I,F,F,X,X,X,X,X,X'
        index=where(long(slitnamearr) eq long(slitname),ct)
        
        isstar=-1
        minwidth=-1
        if ct eq 1 then begin
          pixscale=0.1799
          midpoint=ny/2
          yexpect=midpoint+offsetarr[index[0]]/pixscale
          PRINT,'yexpect:',yexpect
          ;draw white line
          big_img[*,yexpect]=255
        endif else print,'no object found in the slitlist file?!?!?!?'
      endif else print,'no slitlist found for this mask: ',slitlistfile
      
      ;set min/max vals for the image display
      highval=max(big_img)
      lowval=min(big_img)
      

      ;stop
      
      
      extrainfo1=[$
        'CRVAL1',$
        'CDELT1',$
        'CRPIX1',$
        'CTYPE1',$
        'YEXPECT',$
        'MSKNM' ,$
        'FILTNM' ,$
        'SLITNM'  $
        ]
        
      extrainfo2=[$
        string(sxpar(shdr,'CRVAL1')),$
        string(delta),$
        string(refpix),$
        'LINEAR',$
        ssi(yexpect),$
        maskname,$
        filtername,$
        ssi(slitname) $
        ]
        
      ;comments
      extrainfo3=[$
        ' ',$
        ' ',$
        ' ',$
        ' ',$
        ' expected y value of primary obj',$
        ' The mask name', $
        ' The filter name', $
        ' The slit name' $
        ]
        
        
      bmep_display_image,big_img,sciimg,var_img,highval,lowval,slitname,filtername,wavel,savepath,revisevar=0,$
        extrainfo1=extrainfo1,extrainfo2=extrainfo2,extrainfo3=extrainfo3,savetext=1
        
        
    endif ; if files exist
  endfor ; choicearr
  
  
  theend:
  
  cd,original_dir
  print,'end of best mosfire extraction program'
end
