pro bmep_mosdef_check_zspec,zthresh=zthresh

if ~keyword_set(zthresh) then zthresh=0.01
astrolib
cd,getenv('BMEP_MOSDEF_1D')

readcol,'00_redshift_catalog_slim_bmep.txt',masknames,slitnames,ap_nos,z1,ct1,format='A,A,I,F,I,X,X'
z1arr=[]
z1grismarr=[]
zspecarr=[]
zGRISMarr=[]
zPHOTarr=[]
mskarr=[]
slitarr=[]
for i=0,n_elements(masknames)-1 do begin
  if ap_nos[i] eq 1 then begin
    fn=file_search(masknames[i]+'.*.'+slitnames[i]+'*',count=count)
    if count ge 1 then begin
      flux=readfits(fn[0],hdr,/silent,ex=1)
      zspec=sxpar(hdr,'Z_SPEC')
      zgrism=sxpar(hdr,'Z_GRISM')
      zphot=sxpar(hdr,'Z_PHOT')
      if zspec gt 0 then begin
        zspecarr=[zspecarr,zspec]
        z1arr=[z1arr,z1[i]]
        mskarr=[mskarr,masknames[i]]
        slitarr=[slitarr,slitnames[i]]
        endif
      IF zgrism gt 0 then begin
        zGRISMarr=[zGRISMarr,zgrism]
        z1grismarr=[z1grismarr,z1[i]]
        endif
      endif
    endif  
  endfor
plot,z1arr,zspecarr-z1arr,psym=6,/ynozero,ytitle='z spec - z measured',xtitle='z spec'
;plot,zGRISMarr,z1arr,psym=6,/ynozero
;oplot,findgen(100)
  print,'maskname slitname z_spec z_measured'
  forprint,mskarr+' ',slitarr,zspecarr,z1arr
  print
index=where(abs(zspecarr-z1arr) gt zthresh,ct)
if ct ge 1 then begin
  print,'Discrepant  spectra:'
  print,'maskname slitname z_spec z_measured'
  forprint,mskarr[index]+' ',slitarr[index],zspecarr[index],z1arr[index]
  endif

end