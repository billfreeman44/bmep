
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


