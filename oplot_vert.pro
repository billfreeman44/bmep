pro oplot_vert,xpos,yrange,lstyle,name=name,ynamepos=ynamepos,charsize=charsize,thick=thick,cg=cg
color=-1
if keyword_set(name) then if name eq 'Lya' or name eq 'CIV' or name eq '' then color=255 else color=-1
if n_params() lt 3 then lstyle=2
if not keyword_set(charsize) then charsize=2.0
if not keyword_set(thick) then thick=1.0
if keyword_set(cg) then $
  cgplot,[xpos,xpos],yrange,linestyle=lstyle,color='black',thick=thick,/overplot else $
oplot,[xpos,xpos],yrange,linestyle=lstyle,color=color,thick=thick

defsysv,'!xyouts_key',exists=etest
if etest eq 0 then defsysv,'!xyouts_key',0
shiftamt=abs(yrange[0]-yrange[1])/35.0
;,Orientation=90.0

if keyword_set(cg) then cgText,xpos,ynamepos+shiftamt*!xyouts_key,name,charsize=charsize else begin

  if keyword_set(ynamepos) then xyouts,xpos,ynamepos+shiftamt*!xyouts_key,name,charsize=charsize,charthick=1 $
    else begin 
      if keyword_set(name) then xyouts,xpos,0+shiftamt*!xyouts_key,name,charsize=charsize,charthick=4
      endelse
  endelse
  
if !xyouts_key ge  3 then !xyouts_key=0 else !xyouts_key=!xyouts_key+1

end