pro bmep_redshift_analyze,mosdef=mosdef

if ~ keyword_set(mosdef) then savepath=getenv('BMEP_MOSFIRE_1D') $
  else savepath=getenv('BMEP_MOSDEF_1D')
  if savepath eq '' then savepath='~/mosfire/output/idl_output/2D/1d_extracted/'
  
cd,savepath,current=original_dir


readcol,savepath+'00_redshift_catalog_bmep.txt',$
  masknames, filternames, slitnames, apnums,$
  zarr, zerrarr, linenames, restwaves, obswaves,$
  format='A,A,A,A,F,F,A,F,F'

zroundarr=round(zarr*100)
print,n_elements(masknames),' elements before duplicate removal'
i=0
while i lt n_elements(masknames)-1 do begin
  index=where(masknames eq masknames[i] and $
    filternames eq filternames[i] and $
    slitnames eq  slitnames[i] and $
    apnums eq apnums[i] and $
    zroundarr eq zroundarr[i] and $
    linenames eq linenames[i],ct)
    if ct gt 1 then remove,index[1:*],masknames, filternames, slitnames, apnums,$
  zarr, zerrarr, linenames, restwaves, obswaves, zroundarr
  i++
  endwhile 

print,n_elements(masknames),' elements after duplicate removal'

forprint,masknames, filternames, slitnames, apnums,$
  zarr, zerrarr, linenames, restwaves, obswaves,$
  textout=savepath+'00_redshift_catalog_bmep.txt',$
  comment="# maskname filter slit ap_no z zerr linename restwave obswave",$
  format='(A20,A4,A14,A4,F10.6,F13.8,A12,F11.3,F11.3)'

openw,lun,savepath+'00_redshift_catalog_slim_bmep.txt',/get_lun

print,'mask slit apno redshift n_lines'
masknames_nodup=masknames[rem_dup(masknames)]
for i=0,n_elements(masknames_nodup)-1 do begin
  index=where(masknames eq masknames_nodup[i],ct)
  zarr_small=zarr[index]
  zerrarr_small=zerrarr[index]
  slitnames_small=slitnames[index]
  zroundarr_small=zroundarr[index]
  apnums_small=apnums[index]
  slitnames_nodup=slitnames_small[rem_dup(slitnames_small)]
  for j=0,n_elements(slitnames_nodup)-1 do begin ; loop through slits
    for jj=1,7 do begin ;loop through possible aperatures
      ;index is everything in this slit. and is the same aperture number.
      index=where(slitnames_small eq slitnames_nodup[j] and apnums_small eq jj,ct)
      if ct gt 0 then begin
        official_redshifts=zroundarr_small[index[rem_dup(zroundarr_small[index])]]
        official_redshifts=official_redshifts[bsort(official_redshifts)]
        ;loop through and delete anything difference of one
        remove_index=[]
        for k=0,n_elements(official_redshifts)-2 do begin
          if official_redshifts[k] EQ official_redshifts[k+1]-1 then remove_index=[remove_index,k]
          endfor
        if n_elements(remove_index) ne 0 then begin
          print,'REMOVING CLOSE INDICES'
          print,remove_index
          print,official_redshifts
          remove,remove_index,official_redshifts
          print,official_redshifts
          print
          endif
        ca=[]
        za=[]
        for k=0,n_elements(official_redshifts)-1 do begin
          ind=where(abs(zroundarr_small[index]-official_redshifts[k]) le 1,ct)
          weights = zerrarr_small[index[ind]]
          weights[where(weights eq 0,/null)]=0.01
          weights=1.0/weights
          print,masknames_nodup[i],slitnames_nodup[j]+' '+ssi(jj),total(zarr_small[index[ind]]*weights)/total(weights),ct, $
            format='(A10,A13,F8.4,I3)'
          ca=[ca,ct]
          za=[za,total(zarr_small[index[ind]]*weights)/total(weights)]
          endfor; k
        ind=reverse(bsort(ca))
        ca=ca[ind]
        za=za[ind]
        forprint,ca,za
        if n_elements(ca) gt 1 then $
          printf,lun,masknames_nodup[i],slitnames_nodup[j]+' '+ssi(jj),za[0],ca[0],za[1],ca[1],$
            format='(A10,A13,F9.4,I3,F9.4,I3)' $
            else printf,lun,masknames_nodup[i],slitnames_nodup[j]+' '+ssi(jj),za[0],ca[0],-1.0,0,$
            format='(A10,A13,F9.4,I3,F9.4,I3)'
        endif; ct - same aperture, same slit
      endfor; jj
      print
    endfor; j
  endfor; i
close,lun
free_lun,lun

cd,original_dir
print,'program over'
end
