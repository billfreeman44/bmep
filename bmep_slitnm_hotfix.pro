pro bmep_SLITNM_hotfix
key='SLITNM'
cd,getenv('BMEP_MOSDEF_1D')
;cd,'~/mosdef/1d_spec_new'
fn=file_search('*.1d.fits',count=count)
for i=0,n_elements(fn)-1 do begin
  for j=0,6 do begin
    flux=readfits(fn[i],hdr,ex=j,/silent)
    keyout=sxpar(hdr,key)
    print,fn[i],' key read in: ',keyout
    if valid_num(keyout) then  keyout=long(keyout)
    print,fn[i],' key after fix: ',keyout
    print
    sxaddpar,hdr,key,keyout,/savecomment
    modfits,fn[i],flux,hdr,exten_no=j
    endfor 
  endfor

end