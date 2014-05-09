
;6544 7428m
pro bmep_mosdef_getinfo
  !except=2 ;see division by zero errors instantly.
  astrolib
  prior_width=0
  savepath=getenv('BMEP_MOSDEF_1D')
  if savepath eq '' then savepath='~/mosfire/output/idl_output/2D/1d_extracted/'
  cd,savepath,current=original_dir
  
  spawn,'ls *.1d.fits > filenames.txt'
  readcol,'filenames.txt',filenames,format='A'
  spawn,'rm filenames.txt'
  maskarr=[]
  filtarr=[]
  slitarr=[]
  widtharr=[]
  yposarr=[]
  yexpectarr=[]
  objnumarr=[]
  isstararr=[]
  wbyhandarr=[]
  yexpectarr=[]
  minwarr=[]
  blindarr=[]
  
  for i=0,n_elements(filenames)-1 do begin
    data=readfits(filenames[i],hdr,exten_no=1,/silent)
    substrings=strsplit(filenames[i],'.',/extract)
    if n_elements(substrings) eq 5 or n_elements(substrings) eq 6 then begin
      maskarr=[maskarr,sxpar(hdr,'MSKNM')]
      filtarr=[filtarr,sxpar(hdr,'FILTNM')]
      slitarr=[slitarr,sxpar(hdr,'SLITNM')]
      widtharr=[widtharr,sxpar(hdr,'WIDTH')]
      yposarr=[yposarr,sxpar(hdr,'YPOS')]
      objnumarr=[objnumarr,sxpar(hdr,'OBJNUM')]
      isstararr=[isstararr,sxpar(hdr,'ISSTAR')]
      wbyhandarr=[wbyhandarr,sxpar(hdr,'WBYHAND')]
      yexpectarr=[yexpectarr,sxpar(hdr,'YEXPECT')]
      minwarr=[minwarr,sxpar(hdr,'MINW')]
      blindarr=[blindarr,sxpar(hdr,'BLIND')]
    endif
  endfor
  
  
  index=bsort(filtarr)
  maskarr=maskarr[index]
  filtarr=filtarr[index]
  slitarr=slitarr[index]
  widtharr=widtharr[index]
  yposarr=yposarr[index]
  objnumarr=objnumarr[index]
  isstararr=isstararr[index]
  wbyhandarr=wbyhandarr[index]
  yexpectarr=yexpectarr[index]
  minwarr=minwarr[index]
  blindarr=blindarr[index]
  
  index=bsort(objnumarr)
  maskarr=maskarr[index]
  filtarr=filtarr[index]
  slitarr=slitarr[index]
  widtharr=widtharr[index]
  yposarr=yposarr[index]
  objnumarr=objnumarr[index]
  isstararr=isstararr[index]
  wbyhandarr=wbyhandarr[index]
  yexpectarr=yexpectarr[index]
  minwarr=minwarr[index]
  blindarr=blindarr[index]
  
  index=bsort(slitarr)
  maskarr=maskarr[index]
  filtarr=filtarr[index]
  slitarr=slitarr[index]
  widtharr=widtharr[index]
  yposarr=yposarr[index]
  objnumarr=objnumarr[index]
  isstararr=isstararr[index]
  wbyhandarr=wbyhandarr[index]
  yexpectarr=yexpectarr[index]
  minwarr=minwarr[index]
  blindarr=blindarr[index]
  
  ;    index=sort(filtarr)
  ;  maskarr=maskarr[index]
  ;  filtarr=filtarr[index]
  ;  slitarr=slitarr[index]
  ;  widtharr=widtharr[index]
  ;  yposarr=yposarr[index]
  ;  objnumarr=objnumarr[index]
  ;  isstararr=isstararr[index]
  ;  wbyhandarr=wbyhandarr[index]
  ;  yexpectarr=yexpectarr[index]
  
  index=bsort(maskarr)
  maskarr=maskarr[index]
  filtarr=filtarr[index]
  slitarr=slitarr[index]
  widtharr=widtharr[index]
  yposarr=yposarr[index]
  objnumarr=objnumarr[index]
  isstararr=isstararr[index]
  wbyhandarr=wbyhandarr[index]
  yexpectarr=yexpectarr[index]
  minwarr=minwarr[index]
  blindarr=blindarr[index]
  
  
  
  
  openw,lun,savepath+'00_extract_info.txt',/get_lun
  printf,lun,'# mask filter slit-#[*] width[*](starwidth) ypos yposdifference[!!!] actual_width[!!!]'
  printf,lun,'# * after object number means this is the STAR'
  printf,lun,'# * after width means the width was edited by HAND'
  printf,lun,'# & after width means blindly extracted'
  printf,lun,'# !!! after yposdiff means this value is more than 4 different than others'
  printf,lun,'# !!! after actual_width means this value is more than 2.0 different than others'
  print,     'mask filter slit-#[*] width[*](starwidth) ypos yposdifference[!!!] actual_width[!!!]'
  print,     '* after object number means this is the STAR'
  print,     '* after width means the width was edited by HAND'
  print,     '& after width means blindly extracted'
  print,     '!!! after yposdiff means this value is more than 4 different than others'
  print,     '!!! after actual_width means this value is more than 2.0 different than others'
  for i=0,n_elements(slitarr)-1 do begin
  
    ;space out different objects
    if i gt 0 and slitarr[i] ne slitarr[i-1] then print
    if i gt 0 and slitarr[i] ne slitarr[i-1] then printf,lun
    
    ;space out masks
    if i gt 0 and maskarr[i] ne maskarr[i-1] then begin
      ;      printf,lun
      ;      print
      printf,lun,'======================================='
      print,'======================================='
      printf,lun
      print
    endif
    
    if isstararr[i]  eq 1 then suffix='*'  else suffix=' '
    if wbyhandarr[i] eq 1 then suffix2='*' else $
      if blindarr[i]   eq 1 then suffix2='&' else suffix2=' '
    if (i gt 0) and $
      (slitarr[i] eq slitarr[i-1]) and $
      (objnumarr[i] eq objnumarr[i-1]) and $
      (isstararr[i] ne 1) and $
      (abs((yposarr[i]-yexpectarr[i]) - (yposarr[i-1]-yexpectarr[i-1])) gt 4.0) $
      then suffix3='!!! ' else suffix3='    '
    if i gt 0 then prior_width=calc_width
    calc_width=( minwarr[i] gt widtharr[i] or minwarr[i] lt 0.0) ? 0.0: sqrt((widtharr[i]^2-minwarr[I]^2)>0.0)
    if (i gt 0) and $
      (slitarr[i] eq slitarr[i-1]) and $
      (objnumarr[i] eq objnumarr[i-1]) and $
      (isstararr[i] ne 1) and $
      (abs(calc_width - prior_width) gt 2.0) $
      then suffix4='!!! ' else suffix4='    '
      
    print,     maskarr[i],filtarr[i],slitarr[i],objnumarr[i],$
      widtharr[i],minwarr[I],yposarr[i],yposarr[i]-yexpectarr[i],calc_width,$
      format='(A7,A2,I8,"-",I1,"'+suffix+'",F7.2,"'+suffix2+'(",F7.2,") ",F7.2,F7.1,"'+suffix3+'",F6.2,"'+suffix4+'")
    printf,lun,maskarr[i],filtarr[i],slitarr[i],objnumarr[i],$
      widtharr[i],minwarr[I],yposarr[i],yposarr[i]-yexpectarr[i],calc_width,$
      format='(A7,A2,I8,"-",I1,"'+suffix+'",F7.2,"'+suffix2+'",F7.2," ",F7.2,F7.1,"'+suffix3+'",F6.2,"'+suffix4+'")
  ;  stop
  endfor
  index=where(wbyhandarr eq 1)
  print,n_elements(index)
  
  close,lun
  free_lun,lun
  
  cd,original_dir
  print,'program over'
end
