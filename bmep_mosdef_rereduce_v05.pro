;calculate new y profile based on where users clicked. outputs: newydata, newydataerr
;
; hdr - 1d fits header
; 
; d2hdr - 2d fits header
; 
; lun - logical unit number for text file
; 
; sciimg- science image
; 
; var_img - variance image
; 
; newydata - output profile
; 
; newydataerr - output profile error.
;
pro bmep_calc_new_yprofile,hdr,D2hdr,lun,sciimg,var_img,newydata,newydataerr
;convert wavel to pixels in new img.

N_BINS=sxpar(hdr,'NBINS')
wavel2d=(sxpar(D2hdr,'CRVAL1')+findgen(n_elements(sciimg[*,0]))*sxpar(D2hdr,'CDELT1'))
left_pixels=[]
right_pixels=[]
cont_arr=[]

left_bins=[]
right_bins=[]
FOR jj=0,N_BINS-1 do begin
  left_bins=[left_bins,sxpar(hdr,'LBINW'+ssi(jj+1))]
  right_bins=[right_bins,sxpar(hdr,'RBINW'+ssi(jj+1))]
endfor

;make arrays of where users clicked and if it was done in cont mode or not.
FOR jj=0,N_BINS-1 do begin
  index=where(abs(wavel2d-left_bins[jj]) eq min(abs(wavel2d-left_bins[jj])))
  left_pixels=[left_pixels,index[0]]
  index=where(abs(wavel2d-right_bins[jj]) eq min(abs(wavel2d-right_bins[jj])))
  right_pixels=[right_pixels,index[0]]
  cont_arr=[CONT_ARR,sxpar(hdr,'CMODE'+ssi(jj+1),count=count)]
  if count eq 0 then cont_arr[-1]=sxpar(hdr,'CONTMODE') ; older version used this instead of CMODE.
  printf,lun,'LBIN'+ssi(jj+1),left_pixels[jj], sxpar(hdr,'LBIN'+ssi(jj+1))
  printf,lun,'RBIN'+ssi(jj+1),right_pixels[jj],sxpar(hdr,'RBIN'+ssi(jj+1))
  printf,lun,'CMODE'+SSI(JJ+1),cont_arr[jj]
endfor

;collapse columns based on wavelength
newydata=replicate(0.d,n_elements(sciimg[0,*]))
newydataerr=replicate(0.d,n_elements(sciimg[0,*]))
FOR jj=0,N_BINS-1 do begin
  x0=left_pixels[jj]
  x=right_pixels[jj]
  for k=0,n_elements(newydata)-1 do begin
    if x0 lt x then index=indgen(abs(x-x0)>1)+x0 else index=indgen(abs(x-x0)>1)+x
    
    if cont_arr[jj] eq 1 then begin
      noiseguess=[]
      ny=n_elements(var_img[0,*])
      for j=0,n_elements(index)-1 do noiseguess=[noiseguess,total(var_img[index[j],20:ny-21])]
      goodpix=where(noiseguess lt median(noiseguess),ct)
      if ct ge 1 then index=index[goodpix]
    endif
    for j=0,n_elements(index)-1 do newydata[k]=newydata[k]+sciimg[index[j],k]
    for j=0,n_elements(index)-1 do newydataerr[k]=newydataerr[k]+var_img[index[j],k]
  endfor;k=0,n_elements(newydata)-1
endfor;jj=0,N_BINS-1


end



function bmep_slinm_from_filename_1d,filename,mask=mask,apno=apno
  substrings=strsplit(filename,'.',/extract)
  if keyword_set(mask) then return,substrings[0]
  if keyword_set(apno) then begin
    if n_elements(substrings) eq 5 then return,1 else return,fix(substrings[3])
    endif
  return,substrings[2]
end

function bmep_getz_filename,tbl,filename
;APERTURE_NO
;MASKNAME
;SLITOBJNAME
;zmosfire
index=where(tbl.maskname eq bmep_slinm_from_filename_1d(filename,/mask) $
  and tbl.slitobjname eq bmep_slinm_from_filename_1d(filename) $
  and ssi(tbl.APERTURE_NO) eq ssi(bmep_slinm_from_filename_1d(filename,/apno)),ct)
  if ct eq 1 then return,tbl[index].z_mosfire
return,-1.0
end

;;;check which files are different between the old 1d and new 1d folders.
pro bmep_rereduce_sanity_check,twodfolder=twodfolder,onedfolder=onedfolder,tblpath=tblpath
  if ~keyword_set(twodfolder) then path_to_output='/Users/bill/mosdef/00_rereduce_2D/' else path_to_output=twodfolder
  if ~keyword_set(twodfolder) then savepath='/Users/bill/mosdef/00_rereduce_1D/' else savepath=onedfolder
  rereduced_folder_name='000_rereduced'
  
  cd,savepath
  filenames = file_search('*.S*.1d.fits')
  FORPRINT,filenames
  
  filenames = file_search('*.1d.fits')
  print,'list contains ',n_elements(filenames),' objects
  
  
  
  cd,rereduced_folder_name
  re_filenames=file_search('*.1d.fits')

  for i=0,n_elements(filenames)-1 do begin
    index=where(re_filenames eq filenames[i],ct)
    if ct eq 0 then begin
      substrings=strsplit(filenames[i],'.',/extract)
      print,substrings[0]+'-'+substrings[2]+'-'+substrings[3]
      endif
    endfor

;  for i=0,n_elements(re_filenames)-1 do begin
;    index=where(filenames eq re_filenames[i],ct)
;    if ct eq 0 then print,re_filenames[i]+' was not found in the old folder'
;    endfor

end


;create a comparison document for bmep spectra
pro bmep_rereduce_compare,tblpath=tblpath

  if ~keyword_set(tblpath) then tblpath='~/mosdef/Measurements/mosdef_0d_latest.fits'
  savepath='/Users/bill/mosdef/00_rereduce_1D/'
  cd,savepath
  filenames = file_search('*.1d.fits')
  tbl=mrdfits(tblpath,1,/silent)
  currentmask=sss(bmep_slinm_from_filename_1d(filenames[0],/mask))
  warr1=[]
  warr2=[]
  yex1=[]
  yex2=[]
  ypos1=[]
  ypos2=[]
  maskarr=[]
  apnoarr=[]
  fnfinal=[]
  rereduced_folder_name='000_rereduced'
  rereduced_folder_name=rereduced_folder_name+'/'

  CGPS_open,rereduced_folder_name+'0_Extracted_line_comparison_'+currentmask+'.ps',/quiet
  
  for i=0,n_elements(filenames)-1 do begin
;  for i=0,300 do begin
    
    if currentmask ne bmep_slinm_from_filename_1d(filenames[i],/mask) then begin
      cgps_close
      currentmask=sss(bmep_slinm_from_filename_1d(filenames[i],/mask))
      CGPS_open,rereduced_folder_name+'0_Extracted_line_comparison_'+currentmask+'.ps',/quiet
      
      endif
    
    
    data=readfits(filenames[i],hdr,exten_no=1,/silent)
    err=readfits(filenames[i],hdr,exten_no=2,/silent)
    nele=n_elements(data)
    wavel=(sxpar(hdr,'CRVAL1')+findgen(n_elements(data))*sxpar(hdr,'CDELT1'))
    ind_sort=bsort(data)
;    cgplot,wavel,data,yr=[data[ind_sort[nele*5/100]],data[ind_sort[nele*95/100]]]
    
    if file_test(rereduced_folder_name+filenames[i]) then begin
      warr1=[warr1,sxpar(hdr,'WIDTH')]
      yex1=[yex1,sxpar(hdr,'YEXPECT')]
      ypos1=[ypos1,sxpar(hdr,'YPOS')]
      data2=readfits(rereduced_folder_name+filenames[i],hdr,exten_no=1,/silent)
      err2=readfits(rereduced_folder_name+filenames[i],hdr,exten_no=2,/silent)
      wavel2=(sxpar(hdr,'CRVAL1')+findgen(n_elements(data))*sxpar(hdr,'CDELT1'))
      warr2=[warr2,sxpar(hdr,'WIDTH')]
      yex2=[yex2,sxpar(hdr,'YEXPECT')]
      ypos2=[ypos2,sxpar(hdr,'YPOS')]
      maskarr=[maskarr,currentmask]
      apnoarr=[apnoarr,bmep_slinm_from_filename_1d(filenames[i],/apno)]
      fnfinal=[fnfinal,filenames[i]]
;      cgplot,wavel2,data2,/overplot,color='red'
    z=bmep_getz_filename(tbl,filenames[i])
    if z gt 0 then begin
      wavel=wavel/(1.0+z)
      wavel2=wavel2/(1.0+z)
;      linenames= ['HA6565','NII6585','OII_dbl','OIII5008','OIII4960','HB4863']
;      linewavels=[6564.614,6585.27,3728.48,5008.239,4960.295,4862.721];vacuum
      linenames= ['HA6565','OIII5008']
      linewavels=[6564.614,5008.239];vacuum
      for j=0,n_elements(linenames)-1 do begin
        delta=20
        woffset=0.1
        index=where(wavel gt linewavels[j]-10 and wavel lt linewavels[j]+10,ct)
        index2=where(wavel2 gt linewavels[j]-10 and wavel2 lt linewavels[j]+10,ct)
        if ct gt 5 then begin
          cgplot,wavel[index],data[index],xr=[linewavels[j]-10,linewavels[j]+10],title=filenames[i] + ' '+linenames[j],psym=-6
          cgplot,wavel2[index2],data2[index2],/overplot,color='red',psym=-6
          
          cgplot,wavel[index],data[index],xr=[linewavels[j]-10,linewavels[j]+10],title=filenames[i] + ' '+linenames[j],$
            err_yhigh=err[index],err_ylow=err[index]
          cgplot,wavel2[index2]+woffset,data2[index2],/overplot,color='red',$
            err_yhigh=err2[index2],err_ylow=err2[index2]
          endif
          
        
        endfor;linenames
      endif;z gt 0
    endif;file test rereduced
    
    
    
  endfor;filenames
  
  
  
cgps_close

CGPS_open,rereduced_folder_name+'0-compare_widths_ypos_all.ps',/quiet
  cgplot,warr1,warr2,psym=6,xtitle='original width',ytitle='new width'
  cgplot,warr2,warr1-warr2,psym=6,ytitle='original width - new width',xtitle='new width'
  cgplot,ypos1,ypos2,psym=6,xtitle='original ypos',ytitle='new ypos'
  cgplot,yex1-ypos1,yex2-ypos2,psym=6,xtitle='yexpect - ypos original',ytitle='yexpect - ypos new'
  cghistoplot,yex2-ypos2

cgps_close

CGPS_open,rereduced_folder_name+'0-compare_widths_ypos_all_prionly.ps',/quiet
ind=where(apnoarr eq 1)
  cgplot,warr1[ind],warr2[ind],psym=6,xtitle='original width',ytitle='new width'
  cgplot,warr2[ind],warr1[ind]-warr2[ind],psym=6,ytitle='original width - new width',xtitle='new width'
  cgplot,ypos1[ind],ypos2[ind],psym=6,xtitle='original ypos',ytitle='new ypos'
  cgplot,yex1[ind]-ypos1[ind],yex2[ind]-ypos2[ind],$
    psym=6,xtitle='yexpect - ypos original',ytitle='yexpect - ypos new'
  cghistoplot,yex2[ind]-ypos2[ind],xtitle='yexpect - ypos new'
  cghistoplot,yex1[ind]-ypos1[ind],xtitle='yexpect - ypos original'
  
inds=bsort(yex2[ind]-ypos2[ind])
;forprint,yex2[ind[inds]]-ypos2[ind[inds]],' '+fnfinal[ind[inds]]
;print,n_elements(yex2),n_elements(filenames)
cgps_close

mre=rem_dup(maskarr)
CGPS_open,rereduced_folder_name+'0-compare_widths_ypos_all2.ps',/quiet
for i=0,n_elements(mre)-1 do begin
  ind=where(maskarr eq maskarr[mre[i]] and apnoarr eq 1)
;  cgplot,warr1[ind],warr2[ind],psym=6,xtitle='original width',ytitle='new width'
;  cgplot,warr1[ind]-warr2[ind],psym=6,ytitle='original width - new width'
;  cgplot,ypos1[ind],ypos2[ind],psym=6,xtitle='original ypos',ytitle='new ypos'
  cgplot,yex1[ind]-ypos1[ind],yex2[ind]-ypos2[ind],psym=6,$
    xtitle='yexpect - ypos original',ytitle='yexpect - ypos new',title=maskarr[mre[i]],$
    xr=[-3,3],yr=[-3,3]
;  cghistoplot,yex2[ind]-ypos2[ind]
  endfor
cgps_close



end





pro bmep_rereduce_fitgauss,ypos,mwidth,newydata,newydataerr,status,chisq,coeff,gauss_sigma,upper,lower,yfit,mom1,mom2,dummy
        
lower=round(ypos-mwidth)
if lower lt 0 then lower = 0
if lower ge n_elements(newydata) then lower = n_elements(newydata)-1
upper=round(ypos+mwidth)
if upper ge n_elements(newydata) then upper = n_elements(newydata)-1
if upper lt 0 then upper = 0
if upper ne lower then begin
  yfit=newydata[lower:upper]
  yfiterr=sqrt(newydataerr[lower:upper])
  nterms=4
  pi =[{fixed:0, limited:[1,1], limits:[0.D,max(yfit)*1.2]},$ ;peak value
    {fixed:0, limited:[1,1], limits:[-3.D,n_elements(newydata)+2]},$ ;peak centroid
    {fixed:0, limited:[0,0], limits:[0.D,0.D]},$ ;sigma
    {fixed:1, limited:[0,0], limits:[0.D,0.D]}];,$ ;linear bkgnd term
  estimates=[ max(yfit),ypos,2.0,0.0]
  xfit=findgen(upper-lower+1)+lower
  dummy=MPFITPEAK(double(xfit),double(yfit),$
    coeff,nterms=nterms,error=yfiterr,sigma=gauss_sigma,/gaussian,$
    estimates=estimates,parinfo=pi,status=status,chisq=chisq, ERRMSG=errmsg)
  mom1=total(xfit*yfit)/total(yfit)
  mom2=2.355*sqrt(abs(total(xfit*xfit*yfit)/total(yfit) - $
    (total(xfit*yfit)/total(yfit))*(total(xfit*yfit)/total(yfit))))
  chisq=chisq/float(n_elements(yfit)-1.0)
  endif else begin
    status=0
    endelse       
end

function get_primary_yexpect,filename,D2hdr,ny,savepath

  substrings=strsplit(filename,'.',/extract)
  maskname=substrings[0]
  filtername=substrings[1]

;calculate where object SHOULD be!
;note 'shdr' was the science header read in earlier.
yexpect=-1
isstar=-1
minwidth=-1
yshift=0.0
  pixscale=0.1799
  midpoint=ny/2.0
  yexpect=midpoint+sxpar(D2hdr,'OFFSET')/pixscale
  ;check if object is a star, If it isn't shift by amount that star is offset by
  if abs(sxpar(D2hdr,'PRIORITY')) eq 1 then isstar=1 else begin
    isstar=0
    readcol,savepath+'00_starinfo.txt',maskstar,$
      filtstar,objstar,yexpect_star,yactual_star,widthstar,sigmastar,/silent,format='A,A,A,F,F,F,F'
    index=where(maskstar eq maskname and filtstar eq filtername,ct)
    if ct ne 0 then begin
      ;print,ct,' number of stars found for ',maskname,' ',filtername
      yshift=avg(yexpect_star[index] - yactual_star[index])
      minwidth=min(widthstar[index])
    endif else for ii=0,13 do print,'WARNING, DO THE STARS BEFORE DOING ANY OBJECTS'
  endelse ;if isstar
  if sxpar(D2hdr,'SLIT') eq 1 then begin 
    ;print,'this object is at the bottom slit. accounting for this by shifting yexpect by -4'
    yexpect=yexpect-4 ; account for the bottom slit.
    endif
  ;print,'slitname, yexpect, midpoint, yshift, minwidth'
  ;print, yexpect, midpoint, yshift, minwidth
  yexpect=yexpect-yshift
  ;PRINT,'new yexpect:',yexpect



return,yexpect

end



pro bmep_mosdef_rereduce_save,hdr,D2hdr,hdr0,$
  newwidth,minw,ypos,newcenter,mom1,mom2,chisq,coeff,gauss_sigma,$
  i,filenames,fopt,fopterr,f,ferr,newydata,newydataerr,rereduced_folder_name,$
  ny,savepath

extrainfo1=[$
  'CRVAL1',$
  'CDELT1',$
  'CRPIX1',$ 
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
  'MAGNITUD',$
  'PRIORITY',$
  'SCALING',$
  $
  'PA',$
  'CATALOG',$
  'FIELD',$
  'VERSION',$
  'Z_PHOT',$
  $
  'Z_GRISM',$
  'Z_SPEC',$
  'SEEING', $
  'SLITWIDT',$
  'C_OFFSET'$
  ]
  
  extrainfo3=[$
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
  ' Magnitude from MAGMA input files (H band)',$
  ' Priority used in MAGMA   ',$
  ' Scaling factor from cts/s to erg/s/cm^2/Angstrom',$
  $
  ' Slit position angle ',$
  ' Catalog version' , $
  ' ', $
  ' Version of the 3D-HST catalogs', $
  ' Photometric redshift from 3D-HST', $
  $
  ' GRISM redshift from 3D-HST  ', $
  ' External spectroscopic redshift', $
  ' FWHM measured by 2D reduction code [arcsec]', $
  ' Slit width [arcsec]', $
  ' Offset corrected by star position [arscec]' $
  ]
  
 
sxaddpar,hdr,'COMMENT',' Exten 6: Light profile error bars' 
sxaddpar,hdr0,'COMMENT',' Exten 6: Light profile error bars' 
for k=0,n_elements(extrainfo1)-1 do begin 
  sxaddpar,hdr,extrainfo1[k],sxpar(D2hdr,extrainfo1[k],COUNT=COUNT),extrainfo3[k]
  if count ne 1 then print,count,extrainfo1[k]
  endfor
sxaddpar,hdr,'TARGNAME',sxpar(hdr,'TARGNAME')
sxaddpar,hdr,'WIDTH',newwidth,/savecomment
sxaddpar,hdr,'MINW',minw,/savecomment

;;get_primary_yexpect,filename,D2hdr,ny,savepath
if bmep_slinm_from_filename_1d(filenames[i],/apno) eq 1 then $
  sxaddpar,hdr,'YEXPECT',ypos,/savecomment $
  else sxaddpar,hdr,'YEXPECT',GET_PRIMARY_YEXPECT(filenames[i],D2hdr,ny,savepath),/savecomment


sxaddpar,hdr,'YPOS',newcenter,/savecomment
sxaddpar,hdr,'MOM1',mom1,/savecomment
sxaddpar,hdr,'MOM2',mom2,/savecomment
sxaddpar,hdr,'GAUSRCHI',chisq,' Reduced Chi sqr of gauss fit to Y profile'
sxaddpar, hdr, 'GWIDTH', coeff[2], ' Sigma of gaussian fit'
sxaddpar, hdr, 'GCENT', coeff[1], ' Central pixel of gaussian fit'
sxaddpar, hdr, 'GAMP', coeff[0], ' Amplitude of gaussian fit'
sxaddpar, hdr, 'GLINEAR', coeff[3], ' Linear continuum of gaussian fit'
sxaddpar, hdr, 'GWIDTH_E', gauss_sigma[2], ' 1 sigma Error in Sigma of gaussian fit'
sxaddpar, hdr, 'GCENT_E', gauss_sigma[1], ' 1 sigma Error in CENTRAL pixel of gaussian fit'
sxaddpar, hdr, 'GAMP_E', gauss_sigma[0],  ' 1 sigma Error in Amplitude of gaussian fit'
sxaddpar, hdr, 'GLINE_E', gauss_sigma[3], ' 1 sigma Error in Linear continuum of gaussian fit'
FOR K=1,sxpar(D2hdr,'N_OBS') do begin
  sxaddpar,hdr,'FRAME'+ssi(k),sxpar(D2hdr,'FRAME'+ssi(k))
  sxaddpar,hdr,'WEIGHT'+ssi(k),sxpar(D2hdr,'WEIGHT'+ssi(k))
  sxaddpar,hdr,'SEEING'+ssi(k),sxpar(D2hdr,'SEEING'+ssi(k))
  sxaddpar,hdr,'OFFSET'+ssi(k),sxpar(D2hdr,'OFFSET'+ssi(k))
  endfor
              
if sxpar(hdr,'isstar') eq 1 then  sxaddpar,hdr,'MINW',-1,/savecomment             

for k=0,n_elements(extrainfo1)-1 do sxaddpar,hdr0,extrainfo1[k],sxpar(D2hdr,extrainfo1[k]),extrainfo3[k]
sxaddpar,hdr0,'WIDTH',newwidth,/savecomment
sxaddpar,hdr0,'YEXPECT',ypos,/savecomment
sxaddpar,hdr0,'YPOS',newcenter,/savecomment


writefits,rereduced_folder_name+filenames[i],'',hdr0
writefits,rereduced_folder_name+filenames[i],fopt,hdr,/append
writefits,rereduced_folder_name+filenames[i],fopterr,hdr,/append
writefits,rereduced_folder_name+filenames[i],f,hdr,/append
writefits,rereduced_folder_name+filenames[i],ferr,hdr,/append

sxaddpar,hdr,'CRVAL1',1
sxaddpar,hdr,'CDELT1',1.0
sxaddpar,hdr,'CRPIX1',1
sxaddpar,hdr,'CTYPE1','LINEAR'
writefits,rereduced_folder_name+filenames[i],newydata,hdr,/append
writefits,rereduced_folder_name+filenames[i],newydataerr,hdr,/append






end




;compile bmep first.
pro bmep_mosdef_rereduce_v05,twodfolder=twodfolder,onedfolder=onedfolder,tblpath=tblpath
;stop
  FORWARD_FUNCTION bmep_blind_hdr, bmep_dir_exist, bmep_fit_sky,bmep_find_p_slide, $
    bmep_find_p, bmep_get_slitname, bmep_make_hdr,bmep_sigma_clip, bmep_percent_cut
  !except=1 ;see division by zero errors instantly.
  astrolib
  starttime=systime(/seconds)
  mwidth=8 ; width to fit a gaussian to...
  !p.multi=[0,1,1]
  ;spawn,'rm reex_*'
  
  ;create list of 1d spec
;  path_to_output=getenv('BMEP_MOSDEF_2D')
;  if path_to_output eq '' then message,'run in terminal window pls'
  if ~keyword_set(twodfolder) then path_to_output='/Users/bill/mosdef/00_rereduce_2D/' else path_to_output=twodfolder
  
  ;get where to save output of extraction program
;  savepath=getenv('BMEP_MOSDEF_1D')
;  cd,savepath,current=original_dir
  if ~keyword_set(twodfolder) then savepath='/Users/bill/mosdef/00_rereduce_1D/' else savepath=onedfolder
  cd,savepath
  ;read in list of 1d spectra
;  filenames = file_search('co3_01*.1d.fits');for updating only one
  filenames = file_search('*.1d.fits')
  print,'list contains ',n_elements(filenames),' objects



   rereduced_folder_name='000_rereduced'
   if ~bmep_dir_exist(rereduced_folder_name) then file_mkdir,rereduced_folder_name
   rereduced_folder_name=rereduced_folder_name+'/'
  
  
  
  
  ;open file to write info
  openw,lun,rereduced_folder_name+'00_extract_profile_region_info.txt',/get_lun
  printf,lun,'File created on:'
  printf,lun,systime()
  printf,lun,'order is: old val, new val, old - new
  printf,lun
  openw,lun2,rereduced_folder_name+'00_width_ypos_diff_info.txt',/get_lun
  printf,lun2,'File created on:'
  printf,lun2,systime()
  printf,lun2
  openw,lun3,rereduced_folder_name+'00_undone_info.txt',/get_lun
  printf,lun3,'File created on:'
  printf,lun3,systime()
  printf,lun3
  
  ;open plot to plot to compare profiles
  data=readfits(filenames[0],hdr,exten_no=1,/silent)
  cgps_open, Filename=rereduced_folder_name+'0_profile_compare_'+sss(sxpar(hdr,'MSKNM'))+'.ps',/quiet
  

  blindlist=[]
  ;loop through 1d spectra
  for i=0,n_elements(filenames)-1 do begin
;  for i=0, 100 do begin
    old_mask_name=sss(sxpar(hdr,'MSKNM'))
    
    ;read in file
    temp=readfits(filenames[i],hdr0,exten_no=0,/silent)
    data=readfits(filenames[i],hdr,exten_no=1,/silent)
    ydata=readfits(filenames[i],exten_no=5,/silent)
    ydataerr=readfits(filenames[i],exten_no=6,/silent)
    
    ;start new .ps file for each mask.
    if sss(sxpar(hdr,'MSKNM')) ne old_mask_name then begin
      cgps_close
      cgps_open, Filename=rereduced_folder_name+'0_profile_compare_'+sss(sxpar(hdr,'MSKNM'))+'.ps',/quiet
    endif

    blindlist=[blindlist,sxpar(hdr,'BLIND')]
    ;check AUTOEX/WBYHAND/CBYHAND/BLIND flags
    if  $ ;sxpar(hdr,'WBYHAND') eq 1 or
      sxpar(hdr,'AUTOEX')  eq 0 or $
      sxpar(hdr,'CBYHAND') eq 1 or $
      sxpar(hdr,'BLIND')   eq 1 then begin
      if ~sxpar(hdr,'BLIND') then printf,lun,filenames[i]+' not done because of fits hdr flags.'
      if ~sxpar(hdr,'BLIND') then printf,lun3,filenames[i],' because of fits hdr flags.'
    endif else begin
    
    
      ;check for two objects too close to each other.
    
      ;read in variables from 1d spectra
      ypos=sxpar(hdr,'YPOS')
      ypos_old=ypos
      width=sxpar(hdr,'WIDTH')
      wavel=(sxpar(hdr,'CRVAL1')+findgen(n_elements(data))*sxpar(hdr,'CDELT1'))
      sci_file_name=sxpar(hdr,'FILNM')
      printf,lun,filenames[i],' '+sci_file_name
      print,filenames[i]
      
      ;read in NEW! 2d image
      if file_test(path_to_output+'/'+sci_file_name) then begin
        sciimg=readfits(path_to_output+'/'+sci_file_name,D2hdr,exten_no=1,/silent)
        noise_img=readfits(path_to_output+'/'+sci_file_name,D2hdrnois,exten_no=4,/silent)
        ;clean image
        bmep_clean_imgs,sciimg,noise_img
        var_img=noise_img*noise_img
        
        ny=n_elements(sciimg[0,*])
        nx=n_elements(sciimg[*,0])
        
        ;calculate new y data
        bmep_calc_new_yprofile,hdr,D2hdr,lun,sciimg,var_img,newydata,newydataerr
        
        ;plot y cross sections
        ;        cgplot,ydata,title=filenames[i]+'--'+ssi(sxpar(hdr,'CONTMODE')),/nodata
        ;        cgplot,ydata,thick=6,/overplot
        
        ;read in star info
        readcol,savepath+'00_starinfo.txt',maskstar,$
         filtstar,objstar,yexpect_star,yactual_star,widthstar,sigmastar,/silent,format='A,A,A,F,F,F,F'
        readcol,savepath+'00_starinfo_old.txt',maskstar_old,$
         filtstar_old,objstar_old,yexpect_star_old,yactual_star_old,widthstar_old,sigmastar_old,/silent,format='A,A,A,F,F,F,F'
         
        ;set minimum width based on current star info
        index=where(maskstar eq sss(sxpar(hdr,'MSKNM')) and filtstar eq sss(sxpar(hdr,'FILTNM')),ct)
        if ct eq 0 then message,'err, no star'
        minw=min(widthstar[index])
        print,'minw ',minw
        
        
        ;tweak ypos based on difference between old star info and new star info.
        index2=where(maskstar_old eq sss(sxpar(hdr,'MSKNM')) and filtstar_old eq sss(sxpar(hdr,'FILTNM')),ct2)
        if ct2 eq 0 then message,'err, no star'
        ypos=ypos + (avg(yactual_star[index]) - avg(yactual_star_old[index2]))
        cgplot,[ypos,ypos],[min(newydata),max(newydata)],/overplot,color='red',thick=5


        ypos=float(ypos)
        
        
        
        ;normalize cause of errors.
        newydataerr=newydataerr/(max(newydata)*max(newydata))
        newydata=newydata/max(newydata)
        
        cgplot,newydata,thick=3,color='red',title=filenames[i]
        cgplot,ydata,/overplot,color='black'
;        cgerrplot,findgen(n_elements(newydata)),newydata-sqrt(newydataerr),newydata+sqrt(newydataerr)
        cgplot,[ypos_old,ypos_old],[min(newydata),max(newydata)],/overplot,color='black',thick=5
        cgplot,[ypos,ypos],[min(newydata),max(newydata)],/overplot,color='red',thick=5
        
        
        

        ;fit new gaussian
        bmep_rereduce_fitgauss,ypos,mwidth,newydata,newydataerr,status,chisq,coeff,gauss_sigma,upper,lower,yfit,mom1,mom2,dummy

          if status eq 1 or status eq 3 then begin
            fwhm_fit=2.0*SQRT(2.0*ALOG(2.0))*coeff[2]
            if coeff[1] gt -3 and coeff[1] lt n_elements(newydata)+3 then begin

              
              if minw ne -1 and fwhm_fit lt minw then fwhm_fit=minw
              

              newcenter=coeff[1]
              newwidth=fwhm_fit
              if sxpar(hdr,'WBYHAND') eq 1 then newwidth = float(width)
              
              printf,lun,'Width: ',width,newwidth,width-newwidth
              printf,lun,'Center: ',ypos,newcenter,ypos-newcenter
              if newwidth ne width or newcenter ne ypos then begin
                printf,lun2,filenames[i]
                printf,lun2,'Width: ',width,newwidth,width-newwidth
                printf,lun2,'Center: ',ypos,newcenter,ypos-newcenter
                printf,lun2
              endif
              
              cgplot,findgen(upper-lower+1)+lower,yfit,color='green',/overplot
              cgplot,findgen(upper-lower)+lower,dummy,color='purple',/overplot
              cgplot,[coeff[1],coeff[1]],minmax(yfit),color='purple',/overplot
              
              
              if abs(newwidth - width) le 1 and abs(newcenter - ypos) le 3 then begin
                if sxpar(hdr,'ISSTAR') EQ 0 THEN BEGIN
                
                CGTEXT,0.16,0.85,'EXTRACTED',/normal
                
                ;extract and save
                bmep_extraction_simple_subpixel,sciimg,var_img,newydata,newcenter,newwidth,coeff[1],coeff[2],f,ferr,fopt,fopterr,p,1
                
                p=p/max(p)
                p=p*max(yfit)
                
;                cgplot,p,color='blue',/overplot
                
                
              ;save the files...
              bmep_mosdef_rereduce_save,hdr,D2hdr,hdr0,$
                newwidth,minw,ypos,newcenter,mom1,mom2,chisq,coeff,gauss_sigma,$
                i,filenames,fopt,fopterr,f,ferr,newydata,newydataerr,rereduced_folder_name,$
                ny,savepath
              
         
  
                ENDIF;ISSTAR
              endif else begin ;new width and center OK
                printf,lun,'Bad new width or center. Too different'
                printf,lun3,filenames[i],'Bad new width or center. Too different';,newwidth,width,newcenter,ypos
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
        ;printf,lun3,filenames[i]+' - ERROR FINDING 2D FILE, '+path_to_output+'/'+sci_file_name
      endelse
    endelse; header flags ok
    ;end loop
    printf,lun
    endofloop:
  endfor;looping through objects
  
  ;close file and clean up
  close,lun
  free_lun,lun
  close,lun2
  free_lun,lun2
  close,lun3
  free_lun,lun3
  cgps_close
  !p.multi=[0,1,1]
  
  ;this step must be run.
;  bmep_mosdef_update_yexpect

;;SANITY CHECK: create list of missing spectra in rereduced folder compared to 1D.
cd,rereduced_folder_name
slitlist=[]
masklist=[]
yposarr=[]
yexarr=[]
re_filenames=file_search('*.1d.fits')
for i=0,n_elements(filenames)-1 do begin
;for i=0,100 do begin
  if blindlist[i] eq 0 then begin
    index=where(re_filenames eq filenames[i],ct)
    if ct eq 0 then begin
;      print,filenames[i],' is undone'
;      temp=readfits(filenames[i],hdr0,exten_no=0,/silent)
      data=readfits(filenames[i],hdr,exten_no=1,/silent)
      masklist=[masklist,bmep_slinm_from_filename_1d(filenames[i],/mask)]
      slitlist=[slitlist,bmep_slinm_from_filename_1d(filenames[i])]
      yposarr=[]
      yexarr=[]
      endif
    endif
  endfor
  x=rem_dup(slitlist)
  index=[]
  for i=0,n_elements(x)-1 do if strpos(slitlist[x[i]],'S') eq -1 then index=[index,i]
  y=bsort(masklist[x[index]])
  forprint,masklist[x[index[y]]],slitlist[x[index[y]]]
  forprint,masklist[x[index[y]]],' '+slitlist[x[index[y]]],textout='00_undone_slits.txt'
  
  for i=0,n_elements(slitlist[x[index[y]]])-1 do begin
    spawn,'rm '+masklist[x[index[y[i]]]]+'*'+slitlist[x[index[y[i]]]]+'*.fits'
    print,'rm '+masklist[x[index[y[i]]]]+'*'+slitlist[x[index[y[i]]]]+'*.fits'
    endfor

;cd,rereduced_folder_name
  bmep_rereduce_compare,tblpath=tblpath
  
  print,'there are ',file_lines(rereduced_folder_name+'00_undone_info.txt')-3,' undone'
  print,'out of ',n_elements(filenames),' files
  print
  
  print,'that took ',round(systime(/seconds)- starttime),' seconds
  print,'or ',round(systime(/seconds)- starttime)/60.0,' min
  print,'done rereduciing'
  
  
  
  
end