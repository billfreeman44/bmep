pro bmep_redshift_analyze

savepath=getenv('BMEP_MOSDEF_1D')
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

print,'mask slit apno redshift n_lines'
masknames_nodup=masknames[rem_dup(masknames)]
for i=0,n_elements(masknames_nodup)-1 do begin
  index=where(masknames eq masknames_nodup[i],ct)
  slitnames_small=slitnames[index]
  zroundarr_small=zroundarr[index]
  apnums_small=apnums[index]
  slitnames_nodup=slitnames[rem_dup(slitnames_small)]
  for j=0,n_elements(slitnames_nodup)-1 do begin
    for jj=1,7 do begin
      ;index is everything in this slit.
      index=where(slitnames_small eq slitnames_nodup[j] and apnums_small eq jj,ct)
      if ct gt 0 then begin
        official_redshifts=zroundarr_small[index[rem_dup(zroundarr_small[index])]]
        for k=0,n_elements(official_redshifts)-1 do begin
          ind=where(zroundarr_small[index] eq official_redshifts[k],ct)
          print,masknames_nodup[i],slitnames_nodup[j]+'-'+ssi(jj),(official_redshifts[k]/100.0),ct, $
            format='(A10,A13,F6.2,I3)
          endfor; k
        endif; ct
      endfor; jj
      print
    endfor; j
  endfor; i


cd,original_dir
print,'program over'
end
