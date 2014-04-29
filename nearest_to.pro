function nearest_to,array,value
  index=where(abs(array-value) eq min(abs(array-value)),ct)
  ;print,index
  return,index[0]
  end