

;this is an obsolete program that 
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

