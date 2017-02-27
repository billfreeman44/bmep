
;6544 7428m
pro bmep_mosdef_getinfo,yposthresh=yposthresh,widththresh=widththresh
  if ~keyword_set(yposthresh) then yposthresh = 3.0
  if ~keyword_set(widththresh) then widththresh = 1.5
  !except=2 ;see division by zero errors instantly.
  astrolib
  prior_width=0
  savepath=getenv('BMEP_MOSDEF_1D')
;  if savepath eq '' then savepath='~/mosfire/output/idl_output/2D/1d_extracted/005backup/'
  if savepath eq '' then savepath='/Users/bill/mosdef/sedona_1d_extracted/'
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
      slitarr=[slitarr,sss(sxpar(hdr,'SLITNM'))]
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
  
  
  
  print,savepath+'00_extract_info.txt'
  openw,lun,savepath+'00_extract_info.txt',/get_lun
  printf,lun,'# mask filter slit-# [S/N] width [W/B/N] (starwidth) ypos yposdifference [OK/BAD] actual_width [OK/BAD]'
  printf,lun,'# S after object number means this is the STAR (G=GALAXY)'
  printf,lun,'# G after object number means this is a normal galaxy'
  printf,lun,'# W after width means the width was edited by HAND'
  printf,lun,'# B after width means blindly extracted'
  printf,lun,'# N after width means normally extracted'
  printf,lun,'# BAD after yposdiff means this value is more than '+ssf(yposthresh)+' different than others'
  printf,lun,'# BAD after actual_width means this value is more than '+ssf(widththresh)+' different than others'
  print      ,'# mask filter slit-# [S/N] width [W/B/N] (starwidth) ypos yposdifference [OK/BAD] actual_width [OK/BAD]'
  print      ,'# S after object number means this is the STAR (G=GALAXY)'
  print      ,'# G after object number means this is a normal galaxy'
  print      ,'# W after width means the width was edited by HAND'
  print      ,'# B after width means blindly extracted'
  print      ,'# N after width means normally extracted'
  print      ,'# BAD after yposdiff means this value is more than '+ssf(yposthresh)+' different than others'
  print      ,'# BAD after actual_width means this value is more than '+ssf(widththresh)+' different than others'
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
    
    if isstararr[i]  eq 1 then suffix=' S '  else suffix=' G '
    if wbyhandarr[i] eq 1 then suffix2=' W ' else $
      if blindarr[i]   eq 1 then suffix2=' B ' else suffix2=' N '
    if (i gt 0) and $
      (slitarr[i] eq slitarr[i-1]) and $
      (objnumarr[i] eq objnumarr[i-1]) and $
      (isstararr[i] ne 1) and $
      abs((yposarr[i]-yexpectarr[i]) - (yposarr[i-1]-yexpectarr[i-1])) gt yposthresh $
      then suffix3=' BAD ' else suffix3=' OK '
      
    if i gt 0 then prior_width=calc_width
    calc_width=( minwarr[i] gt widtharr[i] or minwarr[i] lt 0.0) ? 0.0: sqrt(((widtharr[i]/2.355)^2-(minwarr[I]/2.355)^2)>0.0)
    if (i gt 0) and $
      (slitarr[i] eq slitarr[i-1]) and $
      (objnumarr[i] eq objnumarr[i-1]) and $
      (isstararr[i] ne 1) and $
      (abs(calc_width - prior_width) gt widththresh) $
      then suffix4=' BAD ' else suffix4=' OK '
      
    print,     maskarr[i]+' ',filtarr[i],slitarr[i],objnumarr[i],$
      widtharr[i],minwarr[I],yposarr[i],yposarr[i]-yexpectarr[i],calc_width,$
      format='(A9,A2,A9,"-",I1,"'+suffix+'",F7.2,"'+suffix2+'(",F7.2,") ",F7.2,F7.1,"'+suffix3+'",F6.2,"'+suffix4+'")
    printf,lun,maskarr[i]+' ',filtarr[i],slitarr[i],objnumarr[i],$
      widtharr[i],minwarr[I],yposarr[i],yposarr[i]-yexpectarr[i],calc_width,$
      format='(A9,A2,A9,"-",I1,"'+suffix+'",F7.2,"'+suffix2+'",F7.2," ",F7.2,F7.1,"'+suffix3+'",F6.2,"'+suffix4+'")
  ;  stop
  endfor
  index=where(wbyhandarr eq 1)
  print,n_elements(index)
  
  close,lun
  free_lun,lun
  
  cd,original_dir
  print,'program over'
end
