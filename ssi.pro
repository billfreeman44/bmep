;==============================================================
;FUNCTION:ssi
;PURPOSE: prevent having to type strcompress(string(fix(number)),/remove_all)
;         every time it is needed
;
;INPUTS: number: an integer
;
;OUTPUTS: a string with no leading or trailing blank spaces
;
;==============================================================
function ssi,number
  return,strcompress(string(fix(number)),/remove_all)
end