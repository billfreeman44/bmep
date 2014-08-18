function ang_delta_to_z_delta,obs,z,delta_ang=delta_ang
if ~keyword_set(delta_ang) then delta_ang=5.0
rest=obs/(z+1.0)
return,1.0 + z -obs/(rest-delta_ang) 

end



pro bmep_mosdef_check_zspec,zthresh=zthresh

if ~keyword_set(zthresh) then zthresh=0.01
astrolib
cd,getenv('BMEP_MOSDEF_1D')

readcol,'00_redshift_catalog_slim_bmep.txt',masknames,slitnames,ap_nos,z1,ct1,format='A,A,I,F,I,X,X'
;if file_test('/Users/bill/mosdef/Measurements/mosdef_0d.fits') then zdcat=mrdfits('/Users/bill/mosdef/Measurements/mosdef_0d.fits')
z1arr=[]
z1grismarr=[]
z1photarr=[]
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
      IF zphot gt 0 then begin
        zPHOTarr=[zPHOTarr,zphot]
        z1photarr=[z1photarr,z1[i]]
        endif
      endif
    endif  
  endfor
ps_start,'00_redshift_plots.ps'
cgplot,z1arr,zspecarr,psym=6,/ynozero,xr=[1,4],yr=[1,4],/xs,/ys,xtitle='z spec',ytitle='z measured'
cgplot,findgen(100),/overplot
cgplot,z1arr,zspecarr-z1arr,psym=6,/ynozero,ytitle='z spec - z measured',xtitle='z measured'
cgplot,replicate(0.01,1000000),/overplot
cgplot,replicate(-0.01,1000000),/overplot
cgplot,z1arr,zspecarr-z1arr,psym=6,/ynozero,ytitle='z spec - z measured',xtitle='z measured',yr=[-0.03,0.03],/ys
;cgplot,replicate(0.01,1000000),/overplot
;cgplot,replicate(-0.01,1000000),/overplot
xarr=findgen(1000)/100.0
;cgplot,xarr,ang_delta_to_z_delta(10000,xarr),/overplot,color='orange'
;cgplot,xarr,(-1.0)*ang_delta_to_z_delta(10000,xarr),/overplot,color='orange'
cgplot,xarr,ang_delta_to_z_delta(15000,xarr),/overplot,color='red'
cgplot,xarr,(-1.0)*ang_delta_to_z_delta(15000,xarr),/overplot,color='red'
;cgplot,xarr,ang_delta_to_z_delta(20000,xarr),/overplot,color='purple'
;cgplot,xarr,(-1.0)*ang_delta_to_z_delta(20000,xarr),/overplot,color='purple'


index=where(abs(zspecarr-z1arr) lt 0.1)
cghistoplot,zspecarr[index]-z1arr[index],binsize=0.0035


cgplot,zGRISMarr,z1grismarr,psym=6,/ynozero,xr=[1,4],yr=[1,4],/xs,/ys,xtitle='z grism',ytitle='z measured'
cgplot,findgen(100),/overplot
cgplot,z1grismarr,zGRISMarr-z1grismarr,psym=6,/ynozero,ytitle='z grism - z measured',xtitle='z measured'
index=where(abs(zGRISMarr-z1grismarr) lt 0.4)
cghistoplot,zGRISMarr[index]-z1grismarr[index],binsize=0.01


cgplot,zPHOTarr,z1photarr,psym=6,/ynozero,xr=[1,4],yr=[1,4],/xs,/ys,xtitle='z phot',ytitle='z measured'
cgplot,findgen(100),/overplot
cgplot,z1photarr,zPHOTarr-z1photarr,psym=6,/ynozero,ytitle='z phot - z measured',xtitle='z measured'
index=where(abs(zPHOTarr-z1photarr) lt 0.4)
;y=histogram(zPHOTarr[index]-z1photarr[index],binsize=0.04,locations=x)
;print,y
;print
;print,x
;help,x,y
;cgplot,float(x),float(y),psym=10
cghistoplot,zPHOTarr[index]-z1photarr[index],binsize=0.04





  print,'maskname slitname z_spec z_measured'
  forprint,mskarr+' ',slitarr,zspecarr,z1arr
  print
index=where(abs(zspecarr-z1arr) gt zthresh,ct)
if ct ge 1 then begin
  print,'Discrepant  spectra:'
  print,'maskname slitname z_spec z_measured'
  forprint,mskarr[index]+' ',slitarr[index],zspecarr[index],z1arr[index]
  endif
ps_end
end