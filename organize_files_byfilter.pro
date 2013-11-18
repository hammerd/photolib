PRO ORGANIZE_FILES_BYFILTER
;
;
; Purpose:
;	Put all flts in same directory. Create README file.
;
;=====================================================================================================
filters=['F098M','F105W','F110W','F125W','F126N','F127M','F128N','F130N','F132N','F139M', 'F140W', 'F153M', 'F160W','F164N', 'F167N','F218W','F225W','F275W','F280N','F300X','F336W','F343N','F373N','F390M','F390W','F395N','F410M','F438W','F467M','F469N','F475W','F502N','F547M','F555W','F606W','F814W','F850LP','G102','G141','G280']
nfilters = n_elements(filters)

;save the filter name for each flt image
dd=findfile('*flt.fits',count=nfiles)
tmpfilt = strarr(nfiles)
for j=0, nfiles-1 do begin
        hdr=headfits(dd(j),exten=0)
        tmpfilt(j) = strupcase(strtrim(sxpar(hdr, 'FILTER'),2))
endfor

;make folders for each unique filter  (USED THIS TO INITIALIZE FILTERS VARIABLE)
;uniq_filt = tmpfilt(uniq(tmpfilt))
;for k=0, n_elements(uniq_filt)-1 do spawn, 'mkdir '+uniq_filt(k)
;stop

; move the flts to the appropriate fitler folder
for k=0, nfilters-1 do begin
	gd=where(tmpfilt eq filters(k), ngd)
	for l=0, ngd-1 do spawn, 'mv '+dd(gd(l))+' '+filters(k)+'/'	;ALL FILTERS
endfor


END

