pro bmep_view_spec_comments 
cd,getenv('BMEP_MOSDEF_1D')
fn=file_search('*1d.fits',count=count)
for i=0,n_elements(fn)-1 do begin
  flux=readfits(fn[i],hdr,ex=1,/silent)
  ucmean=sxpar(hdr,'UCMEAN',count=countucmean)
  ucom=sxpar(hdr,'UCOMMENT',count=countucom)
  if sss(ucmean) ne 'nothingunusual' and countucmean ge 1 then print,fn[i],' ',ucmean
  if sss(ucom) ne 'NoComment' and countucom ge 1 then print,fn[i],' * ',ucom
  endfor
end