;+
; NAME: bisector_los_velocity
;
; PURPOSE:   
;   Compute line-of-sight (LOS) velicity from bisectors at 10 fixed intensity levels
;
; CALLING SEQUENCE:
;   losv = bisector_los_velocity (cube, wavelengths)
;
; + INPUTS:
;   cube:              1D or 3D cube, consisting intensity values of profile(s).
;                      In case of 3D, the cube should be in the form of [x, y, intensity] 
;   wavelengths:       1D array (same size as in the 1D intensity cube) consisting the wavelength of spectral position in Å
;               -- should spatially be the same size as the datacube
;
; + OPTIONAL KEYWORDS/OUTPUT:
;   bisector:         returns the bisector cube
;   loud:             Illustrates the bisectors
;
; + OUTPUTS:
;   losv:             LOS velocity cube (in km/s) corresponding to the 10 bisector (intensity) levels
;                     NOTE: indices 0 to 9 are from the deepest point of the line (i.e., line core) towards the continuum, respectively.
;
; + CREDITS:
;   Author: Shahin Jafarzadeh, March 2015. 
;   Note: Basesd on a bisector mesearment routine as part of LP (La Palma) librrary, developed at the University of Oslo, Norway.
;-

pro bisector_calculate, sp, bisector=bisector, fitp=fitp, loud=loud

; returns bisector of absorption line (a line passing half-way of the 
; connecting points with same intensities); that is a measure of, e.g.,
; the line asymmetry.
; 
; sp      1D array with normalized spectrum
; fitp    number of points to use for line center fit (should be odd number)

 sp=float(sp)/max(sp)  ; nothing should be > continuum
 if keyword_set(fitp) eq 0 then dp=2 else dp=fitp/2
 if keyword_set(loud) eq 0 then loud=0

 bisector=fltarr(2,10)
 xx=findgen(n_elements(sp))

; parabolic fit to line center
 mi=min(sp, mp)
 if mp lt dp then mp=dp  ; should this profile be flagged as weird?
 if mp gt n_elements(sp)-1-dp then mp=n_elements(sp)-1-dp
 co=poly_fit(xx[mp-dp:mp+dp], sp[mp-dp:mp+dp], 2)
 if co[2] gt 0. then begin  ; make sure we do not continue on some strange profile
     ce=(-1.)*co[1]/2/co[2]     ; center
     ldep=co[0] + co[1] * ce + co[2]*ce^2
     bisector[0,0]=ldep         ; line depth at line center
     bisector[1,0]=ce           ; line center position
     p=where(xx eq ce, count)   ; in case line center is found exactly at a sampling point
     if count lt 1 then begin  
         xxb=[xx, ce] & yyb=[sp, ldep]
         ii=sort(xxb) & xxb=xxb[ii] & yyb=yyb[ii]
     endif else begin
         xxb=xx & yyb=sp
     endelse 
     i1=where(xxb eq ce) & i0=0 & i2=n_elements(yyb)-1 ; index range to split profile in 2
; renormalize profile to [0, max(sp)]
     maxsp=max(yyb)
     yyb=yyb-min(yyb) & yyb=yyb/max(yyb) * maxsp
     if (i1 gt 0) and (i1 ne i2) then begin                 ; line center inside sampling range
         m1=floor(min([max(yyb[i0:i1]),max(yyb[i1:i2])]) * 10.) ; last bisector point
         if m1 gt 0 then begin
             if m1 ge 10 then m1=10 ; max bisector at .9
             dd=(findgen(m1)+1.)/10.
             ll0=interpol(xxb[i0:i1], yyb[i0:i1], dd) ; blue part profile
             ll1=interpol(xxb[i1:i2], yyb[i1:i2], dd) ; red part profile
             bisector[0,1:m1]=dd
             bisector[1,1:m1]=(ll1-ll0)/2. + ll0
         endif 
     endif 

     if loud then begin
         window,2
         !p.multi=[0,2,1]
         plot, sp, psym=-1
         oplot, [1,1]*ce, [!y.crange[0], !y.crange[1]], li=1
         oplot, bisector[1,1:m1], bisector[0,0]+bisector[0,1:m1]*(1.-bisector[0,0]), psym=-1
         xxf=findgen((dp*2+1)*10)/10. + mp-dp     
         oplot, xxf, co[0]+co[1]*xxf+co[2]*xxf^2, li=1
         
         plot, xxb, yyb
         oplot, bisector[1,1:m1], bisector[0,1:m1], psym=-1
         oplot, [1,1]*ce, [!y.crange[0], !y.crange[1]], li=1
         !p.multi=0
     endif
 endif 
 return

end

function bisector_los_velocity, cube, wavelengths, bisectors=bisectors, loud=loud, wait=wait
; returns bisector velocitys of given profile(s)

if n_elements(wait) eq 0 then wait = 0
if n_elements(loud) eq 0 then loud = 0

sizecube = size(cube)
if sizecube[0] ne 3 then begin
    if sizecube[0] eq 1 then begin
        blablacube = fltarr(1,1,sizecube[1])
        blablacube[0,0,*] = cube
        cube = blablacube
    endif else stop
endif

sizecube=size(cube)
nx=sizecube[1]
ny=sizecube[2]
ni=sizecube[3] ; number of wavelength positions ==> intensiy values

nl = n_elements(wavelengths)
if nl ne ni then begin
    print, ' '
    print, ' >> number of points in the profile and the number of wavelength positions must be equal! '
    print, ' '
    stop
endif

losVcube = fltarr(nx, ny, 10)
bisectCube=fltarr(nx, ny, 10)

iProfile = cube                           
iProfilem = fltarr(nl)     
for x=0L, nx-1 do for y=0L, ny-1 do iProfilem = iProfilem+iProfile[x,y,*]
iProfilem = iProfilem/max(iProfilem) ; mean Stokes I profile, averaged over the entire FOV

print, '... calculating bisectors'
for x=0L,nx-1 do begin
    for y=0L, ny-1 do begin
        sptemp=reform(iProfile[x,y,*])
        bisector_calculate, sptemp, bisector=bisectoro, loud=loud
        wait, wait
        bisectCube[x,y,0:9]=bisectoro[1,0:9]
    endfor
endfor

; --->>>> indices 0 to 9 are from the deepest point of the line towards the continuum, respectively

bisecim = bisectCube
inanb = where(~finite(bisecim), /null)
bisecim(inanb) = 0. ; pixels with NAN value to zero

print, '... computing LOS velocities'
ref = mean(bisecim[*,*,0:9]) ; mean bisector averaged over the entire FOV, indicating the reference (rest) position
refid = floor(ref)
specsamp = ABS(wavelengths(refid+1)-wavelengths(refid)) ; spectral sampling -  Å
waveref = float((ref-refid)*specsamp)+wavelengths(refid) ; reference (rest) wavelength position - Å
print, '>>> Reference (rest) wavelength: '+strtrim(waveref,2)+' Å'

for gg=0L, 9 do losVcube[*,*,gg] = ((specsamp*(reform(bisecim[*,*,gg])-ref))/waveref)*299792.458

bisectors = bisecim ; optional output

return, reform(losVcube)

print
print, ' ...... returning LOS-velocity cube'
print

end