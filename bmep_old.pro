pro bmep_old,path_to_dropbox=path_to_dropbox,path_to_output=path_to_output,gzending=gzending
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  !except=2

  
  astrolib
  if ~keyword_set(gzending) then gzending=0
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
  filenames = file_search('',/test_directory)
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
;  if gzending eq 0 then spawn,'ls '+maskname+'_?_*eps.fits > xx__spec_list.txt' $
;  else spawn,'ls '+maskname+'_?_*eps.fits.gz > xx__spec_list.txt'
;  readcol,'xx__spec_list.txt',filenames,format='A',/silent
;  spawn,'rm xx__spec_list.txt
  
  filenames = file_search('*eps*')
  if n_elements(filenames) eq 0 then message,'no fits files found... probably bad fliter choice'
  
  ;get the names of the slits
  slitnames=[]
  for i=0,n_elements(filenames)-1 do $
    slitnames=[slitnames,bmep_get_slitname(filenames[i],maskname,/eps,gzending=gzending)]
  slitnames_nodup=slitnames[rem_dup(slitnames)]

  
  
  ;get which slit to do.
  forprint,indgen(n_elements(slitnames_nodup)),replicate('. ',n_elements(slitnames_nodup)),slitnames_nodup
  print,n_elements(slitnames_nodup),' all files'
  print,'which slit?'
  print,'The blank slit is the full 2d image... dont choose this one'
  print
  choice=0
  read,choice
  if choice lt 0 then goto,theend
  if choice ge n_elements(slitnames_nodup)+1 then goto,theend
  
  if choice eq n_elements(slitnames_nodup) then choicearr=indgen(n_elements(filenames)) $
    else choicearr=where(slitnames eq slitnames_nodup[choice])
  
  for i=0,n_elements(choicearr)-1 do begin
    slitname=slitnames[choicearr[i]]
    epsfile=filenames[choicearr[i]]
    ;check if the files exist
    if gzending eq 0 then $
          ivarfile=strmid(epsfile,0,strlen(epsfile)-strlen('_eps.fits'))+'_ivar.fits' $
      else $
      ivarfile=strmid(epsfile,0,strlen(epsfile)-strlen('_eps.fits.gz'))+'_ivar.fits.gz'
    ;check if these files actually exist
    if ~file_test(epsfile) then print,'WARNING!!!! file '+epsfile+' does not exist'
    if ~file_test(ivarfile) then print,'WARNING!!!! file '+ivarfile+' does not exist'
    print,ivarfile
    print,epsfile
    
    if file_test(epsfile) and file_test(ivarfile) and slitname ne '' then begin
      filtername=strmid(epsfile,strlen(maskname)+1,1)
      
      ;read in files
      sciimg=readfits(epsfile,shdr, /SILENT)
      sciimg=double(sciimg)
      ivar_img=readfits(ivarfile, /SILENT)
      ivar_img=double(abs(ivar_img))
      
      ny=n_elements(sciimg[0,*])
      nx=n_elements(sciimg[*,0])
      
      ;calculate variance image
      index=where(ivar_img ne 0.0,ct)
      var_img=replicate(0.0,nx,ny)
      var_img[index]=1.0/ivar_img[index]
      
      ;calc std img
      index=where(ivar_img ne 0.0,ct)
      std_img=replicate(0.0,nx,ny)
      std_img[index]=sqrt(var_img[index])
      
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
      snrimg=sciimg
      snrimg[index]=sciimg[index]/(std_img[index]/2.5)
      snr2sigCut=snrimg
      index=where(snr2sigCut lt 2.0,/null)
      snr2sigCut[index]=0.0
      index=where(snr2sigCut gt 3.0,/null)
      snr2sigCut[index]=3.0
      
      ;mess with these to change how an image is viewed. (scaling)
      botpercent=10.0
      toppercent=90.0
      
;      big_img=findgen(nx,(ny*4))
      big_img=findgen(nx,(ny*2))
      big_img[*,ny*0:ny*1-1]=bytscl(sciimg,top=255,/NAN,$
        min=bmep_percent_cut(sciimg,botpercent),max=bmep_percent_cut(sciimg,toppercent))   ;science img
      big_img[*,ny*1:ny*2-1]=bytscl(snr2sigCut,top=255,/NAN,$
        min=2.0,max=3.0)   ;snr img
;      big_img[*,ny*1:ny*2-1]=bytscl(var_img,top=255,/NAN,$
;        min=bmep_percent_cut(var_img,botpercent),max=bmep_percent_cut(var_img,toppercent))  ;ivar img
;      big_img[*,ny*2:ny*3-1]=bytscl(snrimg,top=255,/NAN,$
;        min=bmep_percent_cut(abs(sciimg),5.0),max=bmep_percent_cut(abs(sciimg),95.0))   ;snr img
;      big_img[*,ny*3:ny*4-1]=bytscl(snr2sigCut,top=255,/NAN,$
;        min=2.0,max=3.0)   ;snr img
;      big_img[*,ny*4:ny*5-1]=bytscl(ivar_img,top=255,/NAN,$
;        min=bmep_percent_cut(ivar_img,botpercent),max=bmep_percent_cut(ivar_img,toppercent))   ;science img
;      big_img[*,ny*5:ny*6-1]=bytscl(std_img,top=255,/NAN,$
;        min=bmep_percent_cut(std_img,botpercent),max=bmep_percent_cut(std_img,toppercent))   ;science img
        
      ;resize big_img if really big in y direction...
      if ny*2 gt 1000 then begin
        big_img=bytscl(sciimg,top=255,/NAN,$
          min=bmep_percent_cut(sciimg,botpercent),max=bmep_percent_cut(sciimg,toppercent))
      endif
      
      ;find yexpect
      slitlistfile='/Users/bill/mosfire/output/idl_output/2D/00mask_info/non_mosdef/masks/'+$
        maskname+'/'+maskname+'_SlitList.txt'
      yexpect=-1
      if file_test(slitlistfile) then begin
        readcol,slitlistfile,slitnamearr,priorityarr,offsetarr,format='X,X,X,X,X,X,X,X,X,I,F,F,X,X,X,X,X,X'
        index=where(sss(slitnamearr) eq sss(slitname),ct)
        
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
        sss(slitname) $
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