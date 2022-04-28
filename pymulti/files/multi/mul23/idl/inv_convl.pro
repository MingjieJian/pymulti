function inv_convl,lambda_air ,debug=debug
;+
;   lambda_vacuum=inv_convl(lambda_air)
;
;-
if(lambda_air gt 2000) then lambda_vacuum=lambda_air/1.00029d0 $
 else lambda_vacuum=lambda_air
loop:
  error=lambda_air-convl(lambda_vacuum)
  lambda_vacuum=lambda_vacuum+error/1.0029d0
  if(keyword_set(debug)) then print,lambda_vacuum,format='(f20.10)'
if(error gt 1.d-4) then goto,loop
return,lambda_vacuum
end
