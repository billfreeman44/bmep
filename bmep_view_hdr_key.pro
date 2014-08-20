pro bmep_view_hdr_key,key,output
cd,getenv('BMEP_MOSDEF_1D')
;cd,'~/mosdef/1d_spec_new'
if N_PARAMS() ge 2 then openw,lun,output,/get_lun
fn=file_search('*.1d.fits',count=count)
for i=0,n_elements(fn)-1 do begin
  flux=readfits(fn[i],hdr,ex=1,/silent)
  keyout=sxpar(hdr,key,count=ckey)
  print,fn[i],' ',keyout
  if N_PARAMS() ge 2 then printf,lun,fn[i],' ',keyout
  endfor
if N_PARAMS() ge 2 then close,lun
if N_PARAMS() ge 2 then free_lun,lun

end