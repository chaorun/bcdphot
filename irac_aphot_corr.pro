;+
; NAME: IRAC_APHOT_CORR
;
;
;
; PURPOSE:     Correct IRAC observed aperture photometry-derived fluxes for the 
;              pixel phase and array location dependent response functions.  ***NOTE:  the pixel phase response
;              is unity (zero variation across a pixel) for channels 3 and 4 (Cryogenic mission only).
;              
;              The correction is of the following form:
;                 CORRECTED_FLUX = Q * Observed_Flux(x,y) / R_al(i,j) * R_pp(xphase,yphase)
;;                
;;                Here the corrected flux is in Jy or proportional units, 
;                 Q is the cal factor converting observed units to Jy (or 
;                 whatever you have chosen), Observed_Flux(x,y) is the measured
;;                aperture flux at centroid position (x,y), R_al(i,j) is 
;                 the Array Location response function (depends on the integer 
;                 pixel i,j), and R_pp(xphase,yphase) is  the Pixel Phase 
;                 response function (depends on pixel phase xphase,yphase). 
;                 Q and the input parameters of R_al and R_pp are hardcoded in this 
;                 program to the values used in deriving the flux conversion factor
;                 in the IRAC cryogenic (S18.24, keyword /CRYO) and warm (S19.2, keyword /WARM)
;                 pipelines.
;                   
;                 You can override these defaults using the PARAMS keyword.  See description of keyword 
;                 PARAMS for more information on these functions.
;
;              Optional keywords /REACH and /HORA allow you to correct data using previously determined cryogenic
;              functions for pixel phase and array location dependent response (Reach et al. 2005; Hora et al. 2008).
;
; CATEGORY:    Photometry, corrections
;
;
;
; CALLING SEQUENCE:   corrected_flux = IRAC_APHOT_CORR ( OBSERVED_FLUX, X, Y, CHANNEL $
;                                               [,PARAMS=PARAMS] [,/PP_ONLY] [,/AL_ONLY] [,/CRYO] [,/WARM]$
;                                               [,/REACH] [,/HORA] $
;                                               [,/SUBARRAY])
;
;
;
; INPUTS:   OBSERVED_FLUX: array of observed aperture fluxes
;                          corresponding to (X,Y)
;                       X: array of X centroids (float), defined such
;                          that the center of the bottom left pixel is
;                          (0,0).  The bottom left corner of the array is
;                          therefore (-0.5,-0.5) and the top right corner is
;                          (255.5,255.5).  Same number of elements as
;                          OBSERVED_FLUX.
;                       Y: array of Y centroids (float). Same number of
;                          elements as OBSERVED_FLUX.
;                 CHANNEL: IRAC channel number (1-4).  Same number of elements
;                          as OBSERVED_FLUX, or a single-element array or scalar
;                          (if only one element, the program assumes all elements of 
;                          OBSERVED_FLUX are for the same channel).
;                          
;                          If CHANNEL is 3 or 4 and none of /CRYO, /REACH, or /HORA is 
;                          set, then the input data will be returned unchanged.
;
;
; OPTIONAL INPUTS: NONE
;
;
;
; KEYWORD PARAMETERS: 
;             PARAMS - a 14 x 4 element matrix, where PARAMS[*,i] gives the
;                        parameter array for the model, for the ith 
;;                       IRAC channel.  ***NOTE: Normally you should not have to change 
;                        the default values of these parameters.
;
;                  Entries for the ith channel of PARAMS
;;                 (i goes from 0 to 3, corresponding to
;;                  channels 1-4) are as follows:
;;
;;                 PARAMS[0:5,i] - A,B,C,D,E,F, the coefficients
;;                                 of the quadratic array
;;                                 location dependence:
;;                 R_al(i,j) = A + B*(i-127) + C*(j-127) +
;;                             D*(i-127)*(j-127) + E*(i-127)^2 +
;;                             F*(j-127)^2
;;                 Here, i = FIX(X) and j = FIX(Y)
;;                -------------------------------------------------- 
;;                         PARAMS[6:12,i] - the coefficients of the
;;                                        double gaussian pixel phase
;;                                        dependence:
;;                 R_pp(xphase,yphase) =
;;                         deltaF_x*exp(-WRAP(xphase-X0,0.5)^2/(2*sigma_x^2))
;;                       + deltaF_y * exp(-WRAP(yphase-Y0,0.5)^2/(2*sigma_y^2))
;;                       + F0
;;                  PARAMS[6,i] = deltaF_x
;;                  PARAMS[7,i] = deltaF_y
;;                  PARAMS[8,i] = X0
;;                  PARAMS[9,i] = Y0
;;                  PARAMS[10,i] = sigma_x
;;                  PARAMS[11,i] = sigma_y
;;                  PARAMS[12,i] = F0_pp
;;    
;;                  Here, xphase and yphase are the observed pixel
;;                  phase (x-ROUND(x) and y-ROUND(y)),
;;                  deltaF_x and deltaF_y are the peak to baseline
;;                  offsets of the x and y gaussians, X0 and Y0 are
;;                  the central x and y pixel phases of the gaussians,
;;                  sigma_x and sigma_y are the gaussian sigmas, and
;;                  F0_pp is the baseline relative flux (relative flux at
;;                  infinity). The function WRAP(z,0.5) wraps the
;;                  pixel phase offset so that it is periodic at +/-
;;                  0.5 pixels (z=+0.5 becomes z=-0.5, continuing to
;;                  0, and vice versa). The center of the wrap (z=0)
;;                  is the phase function peak.
;;                  The pixel phase response is normalized so that the
;;                  function's average, integrated over a pixel,
;;                  equals 1.
;;                -------------------------------------------------- 
;;                  PARAMS[13,i] = the calibration factor converting observed
;                                  units to Jy.  (In most cases, aperture-photometry
;                                  measurements will already be in Jy, so PARAMS[13,*] = 1.0.)
;;
;
;;          /WARM    - Apply the pixel phase and array location corrections appropriate to the post-cryogenic
;;                     (S19.2) pipeline.  This is the default when neither /CRYO, /REACH, or /HORA are set.
;;          /CRYO    - Apply the pixel phase and array location corrections appropriate to the cryogenic 
;;                     (S18.24) pipeline.
;;          /PP_ONLY - Only correct the data using the pixel phase
;;                     function.
;;          /AL_ONLY - Only correct the data using the array location
;;                     dependent function.
;;          /REACH   - Apply the pixel phase and array location corrections of Reach et al. 2005, PASP, 117, 975
;;;         /HORA    - Apply the pixel phase and array location corrections of Hora et al. 2008, PASP, 120, 1233
;           /SUBARRAY- Interpret (X,Y) as relative to the appropriate subarray frame.                 
;
;
; OUTPUTS:  This function returns the corrected flux, which will be either a scalar or 
;           array, depending on the number of elements of OBSERVED_FLUX.
;
;
;
; OPTIONAL OUTPUTS: NONE
;
;
;
; COMMON BLOCKS: NONE
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:  The pixel phase response was derived from aperture photometry performed with a 3 pixel
;                aperture (3-7 pixel background annulus), and so is best used to correct data from the 
;                same aperture. (A larger aperture will give a less pronounced variation in pixel phase 
;                response, in which case this program may overcorrect your data.)
;                
;                This correction will give valid results only when applied to appropriate data.  The default
;                (/WARM) correction is based on S19.2 data from the post-cryogenic mission.  
;                The /CRYO correction is based on S18.18 cryogenic data.  
;                The /REACH and /HORA corrections are relevant to S18.7 and earlier cryogenic data.
;
;
;
; PROCEDURE: 1) Measure X and Y pixel centroid position(s) for stellar image(s) (recommended method: first moment).
;               X and Y are decimal values defined such that the center of the bottom left pixel of the image has (X,Y) = (0.0,0.0). 
;            2) Perform aperture photometry on the observation(s) to derive flux(es) (recommended routine: aper.pro from IDL Astronomy User's Library).  
;               ***NOTE: This correction is optimized for the 3 pixel radius aperture, and 3-7 pixel sky background annulus.
;            3) Correct measured flux(es) using this function. 
;
;
; EXAMPLE:   IDL> corrected_flux = IRAC_APHOT_CORR( OBSERVED_FLUX, X, Y, CHANNEL )
;
;
;
; MODIFICATION HISTORY: 2010 06 28 - Created as APPLY_PIXPHASE_POSDEP         J. Ingalls (SSC)
;                       2010 07 27 - Incorporated new parameter set, fit
;                                    using per-star scale factor and fixed
;                                    position-dependent function
;                                    normalization.           J. Ingalls (SSC)
;                       2010 08 04 - Revised aperture corrections (mode of all
;                                    measurements), renormalizing array location 
;                                    correction to 256x256 pixels, and revised F*K*. 
;                                    Measurement aperture 3,12-20.   J. Ingalls (SSC)
;                       2010 08 05 - Measurement aperture 3, 3-7.    J. Ingalls (SSC)
;                       2010 08 12 - Better photometry with outlier rejection. J. Ingalls (SSC)
;                       2010 09 17 - Photometry derived after fixing flats.  J. Ingalls (SSC)
;                       2010 11 29 - Changed name to IRAC_APHOT_CORR_CRYO.  
;                                    Rewrote header for readability.         J. Ingalls (SSC)
;                       2010 12 02 - Updated /HORA option to include proper normalization
;                                    of array location dependence, and add pixel phase dependence
;                                    Added /REACH option.                    J. Ingalls (SSC)
;                       2011 02 23 - Added /SUBARRAY option to interpret input pixel values as relative to
;                                    each respective subarray FOV.           
;                                    Fixed bug where we get NaN results for CH3 and 4.  J. Ingalls (SSC)
;                       2011 03 28 - Recomputed parameters using the improved color corrections for standard
;                                    stars, and remeasured photometry.  Parameters were confirmed using monte-carlo 
;                                    simulations.                             J. Ingalls (SSC)
;                       2011 03 29 - Recomputed CH3 and 4 parameters after holding fixed at unity the pixel phase 
;                                    function during the fit.  ALSO, fixed an error that set the "A" parameter for the 
;                                    array location correction.               J. Ingalls (SSC)
;                       2011 04 05 - Recomputed all parameters after removing the subarray mode data from the fit,
;                                    which had poorly determined error bars.  J. Ingalls (SSC)
;                       2011 05 12 - Recomputed all parameters after returning the subarray mode data to the fit,
;                                    which had re-computed photometry, using the correct pixel scale, which is different
;                                    in subarray than full array.                         J. Ingalls (SSC)
;                       2015 10 07 - Changed name to IRAC_APHOT_CORR, deprecating IRAC_APHOT_CORR_CRYO.  
;                                    Computed parameters for warm mission.  Added /WARM and /CRYO keywords to 
;                                    distinguish between different regimes.               J. Ingalls (SSC) 
;                       2016 03 02 - Fixed bug that would cause an error when CHANNEL 
;                                    is a vector.  J. Ingalls (SSC)
;                                   
;                                                
;-

FUNCTION aphot_pixphase_posdep_func,x,y,p,pp_only=pp_only,al_only=al_only,reach=reach,hora=hora
;;
;;    Model aperture flux response of an IRAC array as the product of
;;    a quadratic array location dependence and a double gauss pixel
;;    phase dependence.
;;

A = p[0]
B = p[1]
C = p[2]
D = p[3]
E = p[4]
F = p[5]
deltaFx = p[6]
deltaFy = p[7]
x0 = p[8]
y0 = p[9]
sigma_x = p[10]
sigma_y = p[11]
F0 = p[12]   ; IF /REACH or /HORA is set, this parameter doubles as the "A" pixel phase parameter from Hora et al. (2008) Equation (4)
q = p[13]


IF KEYWORD_SET(pp_only) EQ 1 THEN al_model = 1.0 ELSE BEGIN
   ipos = FIX(x)-127.0
   jpos = FIX(y)-127.0
   al_model = (A + B*ipos + C*jpos + D*ipos*jpos + E*ipos^2 + F*jpos^2)
ENDELSE

IF KEYWORD_SET(al_only) EQ 1 THEN pp_model = 1.0 ELSE BEGIN
   xphase = x-ROUND(x)
   yphase = y-ROUND(y)
   IF KEYWORD_SET(REACH) OR KEYWORD_SET(HORA) THEN BEGIN
;;; Pixel phase correction using radial function or Reach et al (2005) or Hora et al. (2008)
      rphase = SQRT(xphase^2 + yphase^2)
      pp_model = 1 + F0 * (1./SQRT(2*!dpi) - rphase)
   ENDIF ELSE BEGIN
;;; Pixel phase correction using double-gauss function in X and Y

      dx_prime = xphase-x0     ;;; Compute wrapped versions of X and Y pixel phases
      iadj_x = where(abs(dx_prime) GT 0.5,nadj_x)
      IF nadj_x NE 0 THEN  dx_prime[iadj_x] = 1.0 - abs(dx_prime[iadj_x])
      dy_prime = yphase-y0
      iadj_y = where(abs(dy_prime) GT 0.5,nadj_y)
      IF nadj_y NE 0 THEN  dy_prime[iadj_y] = 1.0 - abs(dy_prime[iadj_y])

      pp_model = (deltaFx * exp( -dx_prime^2/(2*sigma_x^2) ) + deltaFy * exp( -dy_prime^2/(2*sigma_y^2) ) + F0)
   ENDELSE
ENDELSE

response = al_model * pp_model/q


RETURN,response
END


FUNCTION irac_aphot_corr,flux,x,y,channel,params=params,pp_only=pp_only,al_only=al_only,reach=reach,hora=hora,$
                              subarray=subarray,cryo=cryo,warm=warm

IF N_ELEMENTS(params) EQ 0 THEN BEGIN
   params = MAKE_ARRAY(14,4,value=0.d0)
 
 ;;; 2011 03 28 Revised color corrections, redone photometry,   
 ;;; 2011 03 29 CH 3 and 4 updated to explicitly remove pixel phase function during fit.                  
 ;;; 2011 04 05 No subarray mode data in fit
 ;;; 2011 05 12 Subarray data returned to fit
   IF KEYWORD_SET(CRYO) OR KEYWORD_SET(REACH) OR KEYWORD_SET(HORA) THEN BEGIN
     params[*,0] = [0.98866790,-2.5460463e-05,-4.5413791e-05,-3.7748392e-07,9.6670990e-07,1.1259718e-06,$
       0.018823169,0.030359022,0.091603768,0.0067795815,0.17107575,0.16949466,0.97909886,1.0000000]
     params[*,1] = [0.97713769,0.00016023137,0.00010671203,6.1546421e-07,2.3478283e-06,1.8726664e-06,$
       0.010250904,0.0091393800,0.040266280,0.12475250,0.17673946,0.27301699,0.98964462,1.0000000]
     params[0:5,2] = [0.98318195,-0.00044891858,5.7375573e-05,3.5613363e-07,1.9036209e-06,1.1912425e-06]
     params[0:5,3] = [0.98239745,0.00020132123,-1.9285260e-05,-3.7193490e-07,1.4509036e-06,1.8131923e-06]

     params[[10,11,12,13],2] = 1.0
     params[[10,11,12,13],3] = 1.0
   ENDIF ELSE BEGIN
      WARM = 1
      params[*,0] = [0.98397385,-8.6387206e-05,3.6310403e-05,-5.6284359e-07,1.6023687e-06,1.3574380e-06,$
                    0.039129220, 0.059108698, 0.15899583, 0.037160298, 0.18512016, 0.17005089, 0.95685578,$
                    1.0000000]
      params[*,1] = [0.97494935,0.00015911599,9.5483763e-05,5.5486492e-08,1.2604664e-06, 3.3644277e-06,$
                     0.017579553, 0.0092593634, 0.14839436, 0.11459924, 0.18020032, 0.15649627, 0.98847611,$
                     1.0000000]
     params[0,2] = 1.0
     params[0,3] = 1.0
     params[[10,11,12,13],2] = 1.0
     params[[10,11,12,13],3] = 1.0
      
   ENDELSE
   
 ENDIF
 ;PRINT,params
IF KEYWORD_SET(HORA) OR KEYWORD_SET(REACH) THEN BEGIN
;   PP_ONLY = 0
;   AL_ONLY = 1
   PARAMS = MAKE_ARRAY(14,4,VALUE=0.d0)
   PARAMS[[0,13],*] = 1.0
;; Use F0 for the "A" pixel phase parameter in Hora et al. (2005) Equation (4)
   PARAMS[12,*] = 0.0  ;(setting this to zero means no correction) 
   PARAMS[12,0] = 0.0535
   IF KEYWORD_SET(HORA) THEN PARAMS[12,1] = 0.0309  ;;; only Hora et al (2005) have a correction for Channel 2 
   PARAMS[0:5,0] = [1.0114 ,   -3.536E-6 ,  -6.826E-5 ,  -1.618E-8 ,  1.215E-6  , 1.049E-6]
   PARAMS[0:5,1] = [1.0138 ,    8.401E-5  ,  3.345E-7 ,   1.885E-7  , 1.438E-6  , 1.337E-6 ]
   PARAMS[0:5,2] = [1.0055 ,   -3.870E-4 ,   4.600E-5 ,   1.956E-7 ,  2.078E-6  , 9.970E-7 ]
   PARAMS[0:5,3] = [1.0054 ,    2.332E-4 ,  -8.234E-5  , -1.881E-7  , 6.520E-7 ,  9.415E-7]
;;; Incorporate the appropriate normalization of the array location dependence for Reach (center of array) or Hora (median over entire array)
   IF KEYWORD_SET(REACH) THEN FOR i = 0,3 DO PARAMS[0:5,i] =  PARAMS[0:5,i]/aphot_pixphase_posdep_func(127.0,127.0,PARAMS[*,i],/al_only)
   IF KEYWORD_SET(HORA) THEN BEGIN
      iarray = findgen(256) # MAKE_ARRAY(256,VALUE=1.0)
      jarray = TRANSPOSE(iarray)
      FOR i = 0,3 DO PARAMS[0:5,i] =  PARAMS[0:5,i]/MEDIAN(aphot_pixphase_posdep_func(iarray,jarray,PARAMS[*,i],/al_only))
   ENDIF
ENDIF 
;;;; Interpret input X and Y as relative to appropriate field of view.
IF KEYWORD_SET(SUBARRAY) THEN position_offset = [ [8,216] , [8,216] , [8,8] , [8,8] ] ELSE position_offset = make_array(2,4,value=0)  
   
IF N_ELEMENTS(flux) GT 1 THEN corrected_flux= MAKE_ARRAY(SIZE=SIZE(flux),/DOUBLE) ELSE corrected_flux = 0.d0
IF N_ELEMENTS(channel) EQ 1 THEN BEGIN
   chan_uniq = channel[0]
   IF N_ELEMENTS(flux) GT 1 THEN chn = FIX(MAKE_ARRAY(SIZE=size(flux),value=channel[0]) ) ELSE chn = channel[0]
   nch = 1
ENDIF ELSE BEGIN
   iuniq = UNIQ(channel,SORT(CHANNEL))
   chan_uniq = channel[iuniq]
   nch = N_ELEMENTS(chan_uniq)
   chn = channel
ENDELSE

FOR i = 0,nch-1 DO BEGIN
   Index = WHERE(chn EQ chan_uniq[i],n)
   IF (chan_uniq[i] EQ 3 OR chan_uniq[i] EQ 4) AND KEYWORD_SET(WARM) THEN BEGIN
      PRINT,'IRAC_APHOT_CORR: Channel 3 or 4 data cannot be corrected using /WARM.  Please set /CRYO.  Returning uncorrected fluxes.' 
      corrected_flux[index] = flux[index]
   ENDIF ELSE IF N NE 0 THEN corrected_flux[index] = flux[index]/APHOT_PIXPHASE_POSDEP_FUNC(X[index]+position_offset[0,chan_uniq[i]-1],$
                                          Y[index]+position_offset[1,chan_uniq[i]-1],PARAMS[*,chan_uniq[i]-1],$
                                          pp_only=pp_only,al_only=al_only,reach=reach,hora=hora)
ENDFOR


RETURN,corrected_flux
END
