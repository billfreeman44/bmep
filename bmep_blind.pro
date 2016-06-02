
;create headers for blindly extracted objects
function bmep_blind_hdr,f,extrainfo1,extrainfo2,extrainfo3,yexpect,width,$
    isstar,objnum,min_width,exten=exten,image=image,no_wave_info=no_wave_info
    
  MKHDR, header, f, exten=exten, image=image
  FOR jj=0,n_elements(extrainfo1) -1 do $
    if VALID_NUM(extrainfo2[jj]) then $
    if float(extrainfo2[jj]) eq fix(extrainfo2[jj]) then $
    sxaddpar, Header,extrainfo1[jj],fix(extrainfo2[jj]),extrainfo3[jj] else $
    sxaddpar, Header,extrainfo1[jj],float(extrainfo2[jj]),extrainfo3[jj] else $
    sxaddpar, Header,extrainfo1[jj],STRCOMPRESS(extrainfo2[jj], /REMOVE_ALL),extrainfo3[jj]
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
  sxaddpar, Header, 'SLITLOSS', 0, ' Flux in slit / Total flux assuming 2D Gaussians'
  sxaddpar, Header, 'GWIDTH', -99.0, ' Sigma of gaussian fit'
  sxaddpar, Header, 'GCENT', -99.0, ' Central pixel of gaussian fit'
  sxaddpar, Header, 'GAMP', -99.0, ' Amplitude of gaussian fit'
  sxaddpar, Header, 'GLINEAR', -99.0, ' Linear continuum of gaussian fit'
  sxaddpar, Header, 'GWIDTH_E', -99.0, ' 1sigma Error in Sigma of gaussian fit'
  sxaddpar, Header, 'GCENT_E', -99.0, ' 1sigma Error in CENTRAL pixel of gaussian fit'
  sxaddpar, Header, 'GAMP_E', -99.0,  ' 1sigma Error in Amplitude of gaussian fit'
  sxaddpar, Header, 'GLINE_E', -99.0, ' 1sigma Error in Linear continuum of gaussian fit'
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
  
  sxaddpar, Header, 'COMMENT',' Exten 1: Optimal extraction (not weighted in blind extraction mode)'
  sxaddpar, Header, 'COMMENT',' Exten 2: Optimal extraction error bars'
  sxaddpar, Header, 'COMMENT',' Exten 3: Boxcar extraction'
  sxaddpar, Header, 'COMMENT',' Exten 4: Boxcar extraction error bars'
  sxaddpar, Header, 'COMMENT',' Exten 5: P'
  sxaddpar, Header, 'COMMENT',' Exten 6: P error bars'
  sxaddpar, Header, 'COMMENT',' This is the blind extraction'
  

  
  return, header
  
end


pro bmep_blind_extract,yexpect,width,ny,sciimg,var_img, $
    $ ;OUTPUTS
    f,ferr,fopt,fopterr,p
    
    
    
  f=[]
  ferr=[]
  fopt=[]
  fopterr=[]
  
  bottomint=fix(yexpect-width+1.0)>1
  bottomremainder=1.0 - ( ( (yexpect-width+1.0)-bottomint) < 1.0)
  topint=fix(yexpect+width)<(n_elements(sciimg[0,*])-2)
  topremainder=((yexpect+width)-topint) < 1.0
  
  xarr_small=indgen(topint-bottomint+1)+bottomint
  xarr_small_wider=findgen(topint-bottomint+3)+bottomint-1
  xarr_big=findgen(ny)

  ;STEP2 (steps from horne 1986 extraction paper... step 1 was flatfielding)
  for i=0,n_elements(sciimg[*,0])-1 do begin
  
    ;STEP3  (calculate sky)
    
    ;step4 (extract spectrum and error)
    ;remainder is OUTWARDS from the center...
    botadd= sciimg[i,bottomint-1] * bottomremainder
    botadderr= var_img[i,bottomint-1] * bottomremainder
    topadd= sciimg[i,topint+1] * topremainder
    topadderr= var_img[i,topint+1] * topremainder
    f=[f,total(sciimg[i,bottomint:topint])+topadd+botadd]
    ferr=[ferr,sqrt(total(var_img[i,bottomint:topint])+botadderr+topadderr)]
    
  endfor ; looping through i, the number of columns...
  
  ;optimal extraction based on horne 1989
  ;step 5 (calculate p)
  
  p=replicate(0.0,n_elements(xarr_big))
  p[xarr_small_wider]=replicate(1.0,n_elements(xarr_small_wider))
  p=double(p)
  if total(p) ne 0 then p=p/total(p)

  
  ;set up vars

  ;loop through columns
  for i=0,n_elements(sciimg[*,0])-1 do begin
    ;recalculate fractions to add to top/bottom of extractions.
    botadd= sciimg[i,bottomint-1] * bottomremainder
    botadderr= var_img[i,bottomint-1] * bottomremainder
    topadd= sciimg[i,topint+1] * topremainder
    topadderr= var_img[i,topint+1] * topremainder
    
    ;step 7 (mask cosmic rays)
    ;calculate range over which we care about cosmic rays
    badPixelMask=replicate(1.0,n_elements(xarr_big))
    
    ;swapped steps 6 and 7 to use better flux estimation after removing cosmic spikes.
    ;step 6 (Revise errors) (dont revise for MOSFIRE)
    var_opt=var_img[i,*]
    
    
    ;step 8
    col_data=[botadd,reform(sciimg[i,xarr_small]),topadd]
    var_data=[botadderr,var_opt[xarr_small],topadderr]
    bp_data=badPixelMask[xarr_small_wider]
    p_data=p[xarr_small_wider]
    p_data=p_data/total(p_data)
    
    ;try to not divide by 0.
    index=where(var_data eq 0.0 or var_data eq 99.00,count)
    if count gt 0 then begin
      bp_data[index]=0.0
      var_data[index]=1.0
    endif ;else begin
      
    ;calc optimal flux as in horne
    numerator=total(bp_data*p_data*col_data/var_data)
    denominator=total(bp_data*p_data*p_data/var_data)
    if denominator ne 0 then flux_opt=numerator/denominator else flux_opt=0.0
    ;calc optimal flux error as in horne.
    numerator=total(bp_data*p_data) ;
    denominator=total(bp_data*p_data*p_data/var_data)
    if denominator ne 0 then err_opt=sqrt(abs(numerator/denominator)) else err_opt=0.0
    fopt=[fopt,flux_opt] ;!!!!!!!!!!!!
    fopterr=[fopterr,err_opt]
  endfor ; cols of data! - end of the optimal section of extraction.
   
end

pro bmep_blind_save,savepath,maskname,filtername,slitname,$
    objnum,extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,$
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
  print,'SAVED: '+prefix+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits'
  header=bmep_blind_hdr('',extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,objnum,min_width,/exten)
  writefits,thefilename,'',header
  header=bmep_blind_hdr(fopt,extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,objnum,min_width,/image)
  writefits,thefilename,fopt,header,/append
  header=bmep_blind_hdr(fopterr,extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,objnum,min_width,/image)
  writefits,thefilename,fopterr,header,/append
  header=bmep_blind_hdr(f,extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,objnum,min_width,/image)
  writefits,thefilename,f,header,/append
  header=bmep_blind_hdr(ferr,extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,objnum,min_width,/image)
  writefits,thefilename,ferr,header,/append
  header=bmep_blind_hdr(p,extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,objnum,min_width,/image,/no_wave_info)
  writefits,thefilename,p,header,/append
  header=bmep_blind_hdr(replicate(0.0,n_elements(p)),extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,objnum,min_width,/image,/no_wave_info)
  writefits,thefilename,p,header,/append
endif

;print,maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits'


;if doSave eq 0 then print,'not saved ' else print,'saved'

end



pro bmep_blind,path_to_dropbox=path_to_dropbox,path_to_output=path_to_output
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  starttime=systime(/seconds)
  norepeat=0 ; 1 skips objects already done, 0 redoes objects done blindly.
  !except=2 ;see division by zero errors instantly.
  astrolib
  !p.multi=[0,2,2]
  width=5
  
 ;set output to what is in the envoirnment variable
  x=getenv('BMEP_MOSFIRE_DRP_2D')
  if x ne '' then path_to_output=x
  ;default if no env found
  if ~keyword_set(path_to_output) then path_to_output='~/mosfire/output/idl_output/2D/' ; trailing slash
  
  ;ensure that there is a '/' at the end of the path.
  if strmid(path_to_output,strlen(path_to_output)-1) ne path_sep() then path_to_output=path_to_output+path_sep()
  cd,path_to_output,current=original_dir
  
  ;get where to save output of extraction program
  x=getenv('BMEP_MOSFIRE_DRP_1D')
  if x ne '' then savepath=x else begin
    savepath=path_to_output+'1d_extracted/'
    ;create folder to extract to if it doesn't exist.
    if ~bmep_DIR_EXIST(savepath) then file_mkdir ,'1d_extracted'
    endelse
  ;ensure that there is a '/' at the end of the path.
  if strmid(savepath,strlen(savepath)-1) ne path_sep() then savepath=savepath+path_sep()


  cd,path_to_output,current=original_dir
  
  pwd
  
  
;[path_to_output]/uds_lae5/2014oct3/Y/uds_lae5_Y_n15844_eps.fits
;[path_to_output]/uds_lae5/2014oct29/Y/uds_lae5_Y_n15844_eps.fits
;
;It’s generally like this: [path_to_output]/[maskname]/[date]/[filter]/file.fits
;

;fullfilename=[]
  ;find folders of masks..
  maskfolders = file_search('',/test_directory)
  print,'masks found:'
  forprint,maskfolders
;;;;redundant code because of recursive file_search option. i am a fool.
;  ;loop through mask folder and find dates
;  for i=0,n_elements(maskfolders)-1 do begin
;    cd,maskfolders[i]
;    datefolders = file_search('',/test_directory)
;    print,'dates found:'
;    forprint,datefolders
;    
;    ;loop through date folder and find filters
;    for j=0,n_elements(datefolders)-1 do begin
;      cd,datefolders[j]
;      filterfolders = file_search('',/test_directory)
;      print,'filters found:'
;      forprint,filterfolders
;      
;      ;loop through filter folder and find files
;      for k=0,n_elements(filterfolders)-1 do begin
;        cd,filterfolders[j]
;        filenames = file_search('*_eps.fits',/FULLY_QUALIFY_PATH,count=count)
;        print,ssi(count)+" files found"
;        if count gt 0 then $
;          fullfilename=[fullfilename,filenames]
;          
;        endfor; filterfolders
;        
;      endfor; datefolders
;      
;    endfor;maskfolders

fullfilenames=file_search(x,'*.fits',/full)

filenames=[]
for i=0,n_elements(fullfilenames)-1 do begin
  substrings=strsplit(filenames[i],path_sep(),/extract)
  filenames=[filenames,substrings[-1]]
  endfor

  print,'all 2d files'
  forprint,fullfilename
  print,'all 2d filenames'
  forprint,filenames
  stop
  
  
  
  ;parse folder into different masks!!
  
  
  ;parse names
  masks=[]
  filters=[]
  slitnames=[]
  
  for i=0,n_elements(filenames)-1 do begin
    substrings=strsplit(filenames[i],'_',/extract) ;Note: no '_' in the maskname
    if n_elements(substrings) ge 4 then begin
      masks=[masks,substrings[0:-4]]
      filters=[filters,substrings[-3]]
      slitnames=[slitnames,substrings[-2]]
    endif
  end
  print,'the masks, filters, slitnames found'
  forprint,masks,filters,slitnames
  
  print,'program to stop here'
  stop
  
  ;find 1d extractions of non-primary objects
  print,'creating obj info'
  filenames1d=file_search(savepath+'*.1d.fits' )
  
  openw,lun,savepath+'00_blind_info.txt',/get_lun
  for i=0,n_elements(filenames1d)-1 do begin
    hdr=headfits(filenames1d[i],exten=1)
    objnum=sxpar(hdr,'OBJNUM')
    if sxpar(hdr,'BLIND') ne 1 then begin
      mask=sxpar(hdr,'MSKNM')
      filter=sxpar(hdr,'FILTNM')
      slit=sss(sxpar(hdr,'SLITNM'))
      if valid_num(slit) then slit=ssi(slit)
      width=(sxpar(hdr,'WIDTH'))
      minw=sxpar(hdr,'MINW')
      ypos=sxpar(hdr,'YPOS')
      yexpect=sxpar(hdr,'YEXPECT')
      w_actual_squared=(width*width/(2.355^2) - minw*minw/(2.355^2))>0.0
;      print,slit,' ',w_actual_squared
      printf,lun,mask,filter,slit,objnum,width,$
        w_actual_squared,ypos,yexpect, ypos-yexpect,format='(A10,2x,A4,A10,I3,F9.4,F12.5,F12.5,F9.2,F9.2)' 
      print,mask,filter,slit,objnum,width,$
        w_actual_squared,ypos,yexpect, ypos-yexpect,format='(A10,2x,A4,A10,I3,F9.4,F12.5,F12.5,F9.2,F9.2)' 
    endif
  endfor;n_ele filenames1d
  close,lun
  free_lun,lun
  print,'done creating obj info'
  
;  stop
  
  
  readcol,savepath+'00_blind_info.txt',npmaskarr,npfilterarr,npslitarr,npobjnumarr,$
    npwidtharr,w_actual_sqr_arr,npyposarr,npyexpectarr,npyshiftarr,$
    format="A,A,A,I,I,F,I,F,F"
    
    
  if norepeat eq 0 then PS_Start, Filename=savepath+'00_blind_comparison.ps',/quiet
  
  for i=0,n_elements(slitnames)-1 do begin
    slitname=slitnames[i]
    maskname=masks[i]
    filtername=filters[i]
    filename=maskname+'_'+filtername+'_'+slitname+'_eps.fits'
    noisefilename=maskname+'_'+filtername+'_'+slitname+'_sig.fits'
    
    index=where(npmaskarr eq  maskname and SSS(npslitarr) eq SSS(slitname) and npobjnumarr eq 1,ct)
    if ct ge 1 then w_actual_sqr=avg(w_actual_sqr_arr[index])>0.0 else w_actual_sqr = 0.0

    ;read in files
    sciimg=readfits(filename,shdr, /SILENT,exten_no=1)
    sciimg=double(sciimg)
    
    index=where(finite(sciimg) eq 0,/null)
    sciimg[index]=0.0
    
    ;calculate variance image
    noise_img=readfits(noisefilename, /SILENT,exten_no=1) ;should be a different filename ([...]_sig.fits) with exten_no=1
    
    ;clean image
    index=where(finite(noise_img) eq 0,/null)
    sciimg[index]=0.0
    noise_img[index]=0.0
    noise_img=double(noise_img)
    var_img=noise_img*noise_img
    
    
    ny=n_elements(sciimg[0,*])
    
    ;calculate where object SHOULD be!
    yexpect=-1
    isstar=0
    yshift=0.0
    pixscale=0.1799
    midpoint=ny/2
    yexpect=midpoint+sxpar(shdr,'OFFSET')/pixscale
    ;check if object is a star
    if abs(sxpar(shdr,'PRIORITY')) eq 1 then isstar=1 else  isstar=0
    readcol,savepath+'00_starinfo.txt',maskstar,$
      filtstar,objstar,yexpect_star,yactual_star,widthstar,sigmastar,/silent,format='A,A,A,F,F,F,F'
    index=where(maskstar eq maskname and filtstar eq filtername,ct)
    if ct ne 0 then begin
      ;          print,ct,' number of stars found for ',maskname,' ',filtername
      yshift=avg(yexpect_star[index] - yactual_star[index])
      yexpect=yexpect-yshift
      yexpect=round(yexpect)
      width=(2.355)*sqrt(min(widthstar[index],sub)*min(widthstar[index],sub)/(2.355^2) + w_actual_sqr)
      min_width=min(widthstar[index],sub)
      yposstar=yactual_star[index[sub]]
      starfile=maskname+'.'+filtername+'.'+objstar[index[sub]]+'.1d.fits'
;      p=readfits('1d_extracted/'+starfile,exten_no=5,/silent)
      
      if sxpar(shdr,'SLIT') eq 1 then begin 
        yexpect=yexpect-4 ; account for the bottom slit.
        endif


      
      ;IF THE OBJECT IS EXTRACTED IN OTHER BANDS, THEN USE THEIR WIDTH, NOT THE STAR. (nvm, already fixed)
;      index=where(npmaskarr eq maskname and npslitarr eq slitname AND npobjnumarr eq 1,ct)
;      if ct gt 0 then begin
;        
;        endif;fixing objects 
      
      ;      print,'slitname, yexpect, midpoint, yshift, width'
      print,maskname,' ', filtername,' ', slitname,' ',1, yexpect, midpoint, yshift,min_width, width
            if min_width-0.001 gt width then stop
      
    endif else print,'no object found in the star file?!?!?!?'
    
    
    ;calculate info to add to hdr
    ;its a 2xn array where n is number of things
;    extrainfo1=[$
;      'CRVAL1',$
;      'CDELT1',$
;      'CRPIX1',$
;      'CTYPE1',$
;      'EXPTIME',$
;      'FILNM',$
;      'MSKNM',$
;      'FILTNM',$
;      'SLITNM',$
;      'ISSTAR',$
;      'YEXPECT'$
;      ]
;      
;    extrainfo2=[$
;      string(sxpar(shdr,'CRVAL1')),$
;      string(sxpar(shdr,'CDELT1')),$
;      string(sxpar(shdr,'CRPIX1')),$
;      'LINEAR',$
;      string(sxpar(shdr,'EXPTIME')),$
;      filename,$
;      maskname,$
;      filtername,$
;      slitname, $
;      ssi(isstar), $
;      ssf(yexpect) $
;      ]
      
      
      
       extrainfo1=[$
        'CRVAL1',$
        'CDELT1',$
        'CRPIX1',$
        'CTYPE1',$
        'EXPTIME',$
        $
        'TARGNAME',$
        'MASKNAME',$
        'DATE-OBS',$
        'UT_FIRST',$
        'UT_LAST',$
        $
        'FILTER',$
        'N_OBS',$
        'AIRMASS',$
        'PSCALE',$
        'GAIN',$
        $
        'READNOIS',$
        'PATTERN',$
        'SLIT',$
        'BARS',$
        'RA',$
        $
        'DEC',$
        'OFFSET',$
        'PRIORITY',$
        'SCALING',$
        $
        'PA',$
        'CATALOG',$
        'FILNM',$
        'MSKNM',$
        'FILTNM',$
        $
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
        $
        string(sxpar(shdr,'TARGNAME')),$
        string(sxpar(shdr,'MASKNAME')),$
        string(sxpar(shdr,'DATE-OBS')),$
        string(sxpar(shdr,'UT_FIRST')),$
        string(sxpar(shdr,'UT_LAST')),$
        $
        string(sxpar(shdr,'FILTER')),$
        string(sxpar(shdr,'N_OBS')),$
        string(sxpar(shdr,'AIRMASS')),$
        string(sxpar(shdr,'PSCALE')),$
        string(sxpar(shdr,'GAIN')),$
        $
        string(sxpar(shdr,'READNOIS')),$
        string(sxpar(shdr,'PATTERN')),$
        string(sxpar(shdr,'SLIT')),$
        string(sxpar(shdr,'BARS')),$
        string(sxpar(shdr,'RA')),$
        $
        string(sxpar(shdr,'DEC')),$
        string(sxpar(shdr,'OFFSET')),$
        string(sxpar(shdr,'PRIORITY')),$
        string(sxpar(shdr,'SCALING')),$
        $
        string(sxpar(shdr,'PA')),$
        string(sxpar(shdr,'CATALOG')),$
        filename,$
        maskname,$
        filtername,$
        $
        slitname, $
        ssi(isstar), $
        ssf(min_width), $
        ssf(yexpect)$
        ]
        
      ;comments
      extrainfo3=[$
        ' ',$
        ' ',$
        ' ',$
        ' ',$
        ' Total exposure time (seconds)',$
        $
        ' Name in star list ',$
        ' Name of mask in MAGMA  ',$
        ' Date observed',$
        ' Ut of first obs',$
        ' Ut of first obs',$
        $
        ' Name of the filter',$
        ' Number of included frames',$
        ' Average airmass',$
        ' Pixel scale [arcsec/pix]  ',$
        ' Gain',$
        $
        ' Readnoise',$
        ' Dither pattern',$
        ' Slit Number (Bottom slit is no 1) ',$
        ' Bar numbers (Bottom bar is no 1) ',$
        ' Object Ra (Degrees)',$
        $
        ' Object Dec (Degrees)',$
        ' Offset spectrum wrt center of slit [arscec]',$
        ' Priority used in MAGMA   ',$
        ' Scaling factor from cts/s to erg/s/cm^2/Angstrom',$
        $
        ' Slit position angle ',$
        ' Catalog ',$
        ' name of file',$
        ' name of mask for file naming purposes',$
        ' name of filter for file naming purposes',$
        $
        ' name of slit for file naming purposes', $
        ' Flag if is a star (1 is star, 0 is not)', $
        ' minimum width (-1 default)', $
        ' expected y position (pixels, shifted by star offset)' $
        ]
      
;  FOR K=1,sxpar(shdr,'N_OBS') do begin
;    extrainfo1=[extrainfo1,'FRAME'+ssi(k)]
;    extrainfo2=[extrainfo2,STRING(sxpar(shdr,'FRAME'+ssi(k)))]
;    extrainfo3=[extrainfo3,' ']
;    extrainfo1=[extrainfo1,'WEIGHT'+ssi(k)]
;    extrainfo2=[extrainfo2,STRING(sxpar(shdr,'WEIGHT'+ssi(k)))]
;    extrainfo3=[extrainfo3,' ']
;    extrainfo1=[extrainfo1,'SEEING'+ssi(k)]
;    extrainfo2=[extrainfo2,STRING(sxpar(shdr,'SEEING'+ssi(k)))]
;    extrainfo3=[extrainfo3,' ']
;    extrainfo1=[extrainfo1,'OFFSET'+ssi(k)]
;    extrainfo2=[extrainfo2,STRING(sxpar(shdr,'OFFSET'+ssi(k)))]
;    extrainfo3=[extrainfo3,' Offset in pixels']
;    endfor  

    yexpect=float(yexpect)
    if yexpect ne -1 then begin
    
      bmep_blind_extract,yexpect,width,ny,sciimg,var_img, $
        $ ;OUTPUTS
        f,ferr,fopt,fopterr,p
      objnum=1
      bmep_blind_save,savepath,maskname,filtername,slitname,$
        objnum,extrainfo1,extrainfo2,extrainfo3,yexpect,width,isstar,$
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
      
      index=where(npmaskarr eq maskname and SSS(npslitarr) eq SSS(slitname) and npobjnumarr gt 1,ct)
      if ct gt 0 then begin
        for k=2,6 do begin
          index2=where(npobjnumarr[index] eq k,ct)
          if ct ne 0 then begin ;message,'insanity. probably an object #3 is there with no object #2'
            objnum=k
            npyexpect=round(yexpect+avg(npyshiftarr[index[index2]]))
            npwidth=(2.355)*sqrt(min_width*min_width/(2.355^2) + avg(w_actual_sqr_arr[index[index2]]))
            
            print,maskname,' ', filtername,' ', slitname,' ',k, npyexpect, midpoint, yshift, min_width, npwidth
            if min_width-0.001 gt npwidth then stop
            
            bmep_blind_extract,$
              npyexpect,$ ; yposition
              npwidth,$ ; width
              ny,sciimg,var_img,f,ferr,fopt,fopterr,p
              
            bmep_blind_save,savepath,maskname,filtername,slitname,$
              objnum,extrainfo1,extrainfo2,extrainfo3,npyexpect,npwidth,0,$ ; 0 is the isstar parameter.
              f,ferr,fopt,fopterr,p,min_width
              
            suffix='.'+ssi(objnum)
;            if file_test(savepath+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits')  then begin
;              data=readfits(savepath+maskname+'.'+filtername+'.'+slitname+suffix+'.1d.fits',shdr,exten_no=5,/silent)
;              cgplot,data,title=maskname+'.'+filtername+'.'+slitname+suffix+' y profile',ytitle='P'
;              cgplot,(p/max(p))*max(data),color='red',/overplot
;            endif
            
          endif; ct ne 0
        endfor; k
      endif ; else print,'no np objects' ;ct gt 0
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



