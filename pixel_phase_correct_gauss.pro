;+
; NAME: PIXEL_PHASE_CORRECT_GAUSS
;
;
;
; PURPOSE: Correct for pixel phase the aperture photometry flux(es)
;          OBSERVED_FLUX of point sources observed with Warm IRAC,
;          given a moment centroid at pixel position(s) X,Y.  (The
;          coordinate system here is defined with the center of a
;          pixel having integer values and the edges having
;          half-integer values).  The model of pixel phase response is
;          a "double-gaussian", summation of gaussians in the
;          orthogonal pixel phase directions:
;          
;      RESPONSE(xphase,yphase) = F0 + deltaF_x*exp(-WRAP(xphase-X0,0.5)^2/(2*sigma_x^2)) + 
;                                     deltaF_y*exp(-WRAP(yphase-Y0,0.5)^2/(2*sigma_y^2)) 
;
;          Here, xphase and yphase are the observed pixel phase
;          (X-ROUND(X) and Y-ROUND(Y)), deltaF_x and
;          deltaF_y are the peak to baseline offsets of the x and y
;          gaussians, X0 and Y0 are the central x and y pixel phases
;          of the gaussians, sigma_x and sigma_y are the gaussian
;          sigmas, and F0_pp is the baseline relative flux (relative
;          flux at infinity). The function WRAP(z,0.5) wraps the pixel
;          phase offset so that it is periodic at +/- 0.5 pixels
;          (z=+0.5 becomes z=-0.5, continuing to 0, and vice
;          versa). The center of the wrap (z=0) is the phase function
;          peak.  The response function is normalized such that its
;          integral over a pixel is unity.  In other words, for random
;          pointing, the average correction is 1.
;
;
; CATEGORY: Photometry, corrections
;
;
; CALLING SEQUENCE:  IDL> corrected_flux = PIXEL_PHASE_CORRECT_GAUSS(observed_flux,x,y,channel [,APER]
;                                          [,PARAMS=PARAMS]
;                                          [,/SUBARRAY] [,/HELP])
;
;
;
; INPUTS:  OBSERVED_FLUX  - A scalar or vector containing observed
;                           aperture flux measurements in Jy or
;                           proportional units.
;          X, Y           - Scalars or vectors with the same number of
;                           elements as OBSERVED_FLUX, giving the X
;                           and Y pixel centroid location(s) of the
;                           source for each flux measurement.  For the
;                           correction to be valid, centroids should
;                           be obtained using the moment method.  
;                           Alternately, X and/or Y can be the X and Y
;                           pixel phases of the measurements (the
;                           pixel phase correction does not currently depend on
;                           absolute position on the array, only pixel
;                           phase). The coordinate system is defined
;                           with the center of a pixel having a whole
;                           number pixel value (pixel phase 0.0) and
;                           the edge having half-integer values (pixel
;                           phase +/- 0.5).
;           CHANNEL        -The IRAC channel number for the
;                           observations.  If a scalar value, then all
;                           observations will be assumed to have the
;                           same channel.  If a vector the it should
;                           have the same number of elements as
;                           OBSERVED_FLUX (otherwise only the first
;                           N_ELEMENTS(channel) elements of
;                           OBSERVED_FLUX will be corrected).
;                           Allowable values are 1 and 2.  Can be an
;                           integer, floating point or string
;                           variable.
;
;
; OPTIONAL INPUTS:      APER  -An optional scalar string describing
;                              the aperture used to perform the
;                              centroiding and photometry.  This entry
;                              consists of three integers, separated
;                              by underscores (_): 
;                                (1) the radius of the aperture, in
;                                pixels; 
;                                (2) the inner radius of the sky
;                                annulus, in pixels; and
;                                (3) the outer radius of the sky
;                                annulus.
; 
;                              Currently available values are:
;
;                              APER     RADIUS       ANNULUS 
;                              =============================
;                              2_2_6      2           2 -  6
;                              2_12_20               12 - 20
;                              3_3_7      3           3 -  7
;                              3_12_20               12 - 20
;                              5_5_10     5           5 - 10
;                              4_4_8      4           4 -  8
;                              4_12_20    4          12 - 20
;                              5_12_20               12 - 20
;                              8_12_20    8          12 - 20
;                              9_12_20    9          12 - 20
;                              10_12_20  10          12 - 20   ***DEFAULT***
;                              11_12_20  11          12 - 20
;                              12_12_20  12          12 - 20
;                              15_15_25  15          15 - 25
;   
;                           If APER is not passed, or is not passed with one
;                           of the values above, it will default to '10_12_20'.
;
;
; KEYWORD PARAMETERS:  PARAMS  - An optional structure variable that allows you
;                                to override the default pixel phase
;                                response parameters.  The structure is
;                                initialized/defined using the
;                                following syntax (or similar):
;                  IDL> params = {x0:fltarr(2),y0:fltarr(2),sigma_x:fltarr(2),sigma_y:fltarr(2),$
;                                 deltaF_x:fltarr(2),deltaF_y:fltarr(2),F0:fltarr(2),theta:fltarr(2)}
;
;                                Each tag of PARAMS specifies a
;                                2-element vector, whose 2 elements
;                                give the pixel phase response
;                                function parameters for Channels 1
;                                and 2:
;                      PARAMS.X0:  Peak pixel phase of the X gaussian
;                      PARAMS.Y0:  Peak pixel phase of the Y gaussian
;                      PARAMS.SIGMA_X: Sigma width of the X gaussian  
;                      PARAMS.SIGMA_Y: Sigma width of the Y gaussian  
;                      PARAMS.deltaF_X: Peak-to-baseline offset of
;                                       the X gaussian  
;                      PARAMS.deltaF_Y: Peak-to-baseline offset of
;                                       the Y gaussian  
;                      PARAMS.F0: Baseline response at infinity
;                      PARAMS.THETA: Angle (in degrees) that the X
;                                    Gaussian makes with the X pixel
;                                    axis.  (Y Gaussian is orthogonal
;                                    to X Gaussian).  In the default
;                                    model, the angle is assumed to be
;                                    0.
;                       If you set PARAMS, then APER will be ignored
;
;                       /SUBARRAY - use the parameter set for subarray mode
;
;                       /HELP - print this header
;
;
; OUTPUTS:  This routine returns the corrected flux value(s) in the same
;           units as OBSERVED_FLUX.  The result will have the same
;           number of elements and type as OBSERVED_FLUX.
;
;
;
; OPTIONAL OUTPUTS: None
;
;
;
; COMMON BLOCKS: None
;
;
;
; SIDE EFFECTS:  The response function is normalized such
;      that its integral over a pixel is unity.  In other words, for
;      random pointing, the average result would be no correction.
;
;
;
; RESTRICTIONS:  (1) Inputs OBSERVED_FLUX, X, and Y should have
;                the same number of elements.  If CHANNEL is a vector,
;                then it should also have the same number of elements
;                as OBSERVED_FLUX, X, and Y.  
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY: 2/24/10: Created.  Parameters based on Focus
;                                Check data.            J. Ingalls (SSC)
;                       4/20/10: Introduced PARAMS structure to
;                                override default parameters.  Changed
;                                default parameters to reflect Full
;                                Array measurements, with more robust
;                                normalization.         J. Ingalls (SSC)
;                       5/03/10: Introduced APER parameter to accomodate
;                                variations of pixel phase response with
;                                aperture.  Added /HELP keyword.
;                                                       J. Ingalls (SSC)
;                       5/18/10: Introduced /SUBARRAY parameter to 
;                                correct pixel phase response in
;                                subarray FOV.  
;                                                       J. Ingalls (SSC)
;                       10/27/11: Added 4-pixel radius aperture (with 4,8 and 12,20 background),
;                                 Recomputed parameters for all apertures using
;                                 updated centroiding technique.  Now all apertures have
;                                 valid results.
;                                                       J. Ingalls (SSC)
;
;    Copyright (C) <2011>  <J. Ingalls/SSC>
;    This program is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see <http://www.gnu.org/licenses/>.
;
;-
FUNCTION pixel_phase_double_gauss,xphase,yphase,x0,y0,sigma_x,sigma_y,deltaF_x,deltaF_y,F0,theta

dx = xphase-x0
dy = yphase-y0
IF Theta EQ 0 THEN BEGIN
   dx_prime = dx 
   dy_prime = dy
ENDIF ELSE BEGIN
   dx_prime = dx * cos(theta) + dy * sin(theta)
   dy_prime = -dx * sin(theta) + dy * cos(theta)
ENDELSE
iadj_x = where(abs(dx_prime) GT 0.5,nadj_x)
IF nadj_x NE 0 THEN  dx_prime[iadj_x] = 1.0 - abs(dx_prime[iadj_x])
iadj_y = where(abs(dy_prime) GT 0.5,nadj_y)
IF nadj_y NE 0 THEN  dy_prime[iadj_y] = 1.0 - abs(dy_prime[iadj_y])

relative_flux = deltaF_x * exp( -dx_prime^2/(2*sigma_x^2) ) + deltaF_y * exp( -dy_prime^2/(2*sigma_y^2) ) + F0

RETURN,relative_flux
END

FUNCTION pixel_phase_correct_gauss,observed_flux,x,y,channel,aper,params=params,subarray=subarray,help=help
;;
;; Correct the flux(es) OBSERVED_FLUX of measured objects with moment centroid at
;; pixel position(s) X,Y.  Use the pixel phase response double
;; gaussian fit.  OBSERVED_FLUX,X,Y, and CHANNEL are arrays or scalars with the
;; same number of elements
;;
IF KEYWORD_SET(HELP) THEN doc_library,'pixel_phase_correct_gauss'

AP_ARRAY = ['15_15_25','12_12_20','11_12_20','10_12_20','9_12_20','8_12_20','5_12_20','5_5_10','4_12_20','4_4_8','3_12_20','3_3_7','2_12_20','2_2_6']

IF N_ELEMENTS(APER) EQ 0 THEN aper = '10_12_20'
nap = 0
WHILE nap EQ 0 DO BEGIN
   iap = where(AP_ARRAY EQ APER,nap)
   IF NAP EQ 0 THEN aper = '10_12_20'
ENDWHILE

nflux = N_ELEMENTS(observed_flux)
nx = N_ELEMENTS(x)
ny = N_ELEMENTS(y)
nch = N_ELEMENTS(channel)
IF nch EQ 1 THEN BEGIN
   IF Nflux GT 1 THEN channel = FIX(MAKE_ARRAY(SIZE=size(observed_flux),VALUE=channel)) ELSE channel = FIX(channel)
   nch = N_ELEMENTS(channel)
ENDIF

IF nflux NE nx OR nflux NE ny OR nflux NE nch THEN BEGIN
   print,'PIXEL_PHASE_CORRECT_GAUSS:  Input array sizes don''t match.  Exiting.'
   RETALL
ENDIF

IF N_ELEMENTS(PARAMS) EQ 0 then BEGIN
;;; Parameters constrained such that average over pixel is unity.

   IF KEYWORD_SET(SUBARRAY) EQ 0 THEN BEGIN
   ;;; From run of 10/25/11, adding 4-pixel radius apertures and using new centroiding function
      X0_CH1 = [-0.075,0.192,0.192,0.192,0.191,0.192,0.191,0.191,0.194,0.193,0.190,0.186,0.208,0.210]
      Y0_CH1 = [0.051,0.029,0.029,0.030,0.029,0.030,0.027,0.027,0.023,0.023,0.022,0.021,0.036,0.034]
      SIG_X_CH1 = [0.990,0.175,0.187,0.187,0.189,0.189,0.189,0.190,0.190,0.190,0.189,0.191,0.198,0.197]
      SIG_Y_CH1 = [0.191,0.169,0.164,0.163,0.162,0.161,0.163,0.163,0.163,0.164,0.164,0.164,0.163,0.163]
      DELF_X_CH1 = [0.0396,0.0279,0.0286,0.0288,0.0290,0.0291,0.0299,0.0298,0.0309,0.0308,0.0327,0.0330,0.0351,0.0350]
      DELF_Y_CH1 = [0.0433,0.0451,0.0448,0.0453,0.0453,0.0461,0.0475,0.0479,0.0487,0.0492,0.0544,0.0540,0.0633,0.0634]
      F0_CH1 = [0.941,0.969,0.968,0.968,0.968,0.968,0.967,0.966,0.966,0.965,0.962,0.962,0.957,0.957]
      X0_CH2 = [-0.044,0.089,0.089,0.089,0.088,0.089,0.088,0.091,0.083,0.087,0.070,0.071,0.047,0.051]
      Y0_CH2 = [-0.179,0.100,0.099,0.100,0.098,0.098,0.095,0.091,0.096,0.095,0.098,0.099,0.095,0.104]
      SIG_X_CH2 = [0.014,0.218,0.217,0.215,0.212,0.210,0.209,0.209,0.206,0.207,0.203,0.201,0.100,0.104]
      SIG_Y_CH2 = [0.018,0.227,0.227,0.228,0.224,0.221,0.214,0.209,0.215,0.210,0.216,0.218,0.177,0.191]
      DELF_X_CH2 = [0.1457,0.0172,0.0173,0.0172,0.0170,0.0171,0.0178,0.0176,0.0173,0.0177,0.0175,0.0179,0.0078,0.0083]
      DELF_Y_CH2 = [0.0312,0.0188,0.0189,0.0189,0.0187,0.0185,0.0187,0.0183,0.0188,0.0191,0.0186,0.0192,0.0260,0.0270]
      F0_CH2 = [0.994,0.980,0.980,0.980,0.981,0.981,0.981,0.981,0.981,0.981,0.981,0.981,0.987,0.985]

   ENDIF ELSE BEGIN
   ;;; From run of 10/25/11, adding 4-pixel radius apertures
      X0_CH1 = [0.125,0.125,0.129,0.129,0.130,0.130,0.130,0.129,0.131,0.130,0.128,0.126,0.130,0.131]
      Y0_CH1 = [0.039,0.036,0.033,0.033,0.033,0.033,0.032,0.032,0.030,0.032,0.027,0.027,0.039,0.037]
      SIG_X_CH1 = [0.194,0.187,0.184,0.184,0.184,0.184,0.183,0.185,0.183,0.184,0.183,0.182,0.177,0.176]
      SIG_Y_CH1 = [0.192,0.194,0.191,0.192,0.192,0.193,0.194,0.196,0.194,0.197,0.195,0.197,0.191,0.191]
      DELF_X_CH1 = [0.0358,0.0362,0.0359,0.0360,0.0363,0.0366,0.0384,0.0388,0.0396,0.0394,0.0416,0.0419,0.0425,0.0430]
      DELF_Y_CH1 = [0.0434,0.0438,0.0434,0.0437,0.0438,0.0447,0.0468,0.0472,0.0476,0.0488,0.0520,0.0521,0.0596,0.0596]
      F0_CH1 = [0.962,0.962,0.963,0.963,0.962,0.962,0.960,0.959,0.959,0.958,0.956,0.956,0.953,0.953]
      X0_CH2 = [0.071,0.074,0.070,0.067,0.065,0.064,0.052,0.051,0.048,0.049,0.035,0.037,0.009,0.011]
      Y0_CH2 = [0.106,0.105,0.107,0.113,0.111,0.110,0.099,0.098,0.089,0.088,0.091,0.095,0.131,0.133]
      SIG_X_CH2 = [0.175,0.175,0.172,0.164,0.167,0.168,0.173,0.177,0.187,0.188,0.196,0.198,0.171,0.174]
      SIG_Y_CH2 = [0.265,0.291,0.284,0.267,0.263,0.257,0.212,0.207,0.192,0.180,0.179,0.174,0.177,0.179]
      DELF_X_CH2 = [0.0167,0.0182,0.0184,0.0185,0.0182,0.0182,0.0185,0.0187,0.0193,0.0193,0.0215,0.0220,0.0205,0.0214]
      DELF_Y_CH2 = [0.0169,0.0201,0.0195,0.0189,0.0192,0.0187,0.0167,0.0167,0.0162,0.0161,0.0162,0.0169,0.0223,0.0236]
      F0_CH2 = [0.982,0.979,0.979,0.981,0.980,0.981,0.983,0.983,0.983,0.984,0.982,0.982,0.981,0.980]

   ENDELSE

;;;              CH1         CH2
   x0 = [x0_ch1[iap[0]],x0_ch2[iap[0]]]
   y0 = [y0_ch1[iap[0]],y0_ch2[iap[0]]]
   sigma_x = [sig_x_ch1[iap[0]],sig_x_ch2[iap[0]]]
   sigma_y = [sig_y_ch1[iap[0]],sig_y_ch2[iap[0]]]
   deltaF_x = [delf_x_ch1[iap[0]],delf_x_ch2[iap[0]]]
   deltaF_y = [delf_y_ch1[iap[0]],delf_y_ch2[iap[0]]]
   F0 = [F0_ch1[iap[0]],F0_ch2[iap[0]]]
   Theta = [0.,0.]

ENDIF ELSE BEGIN
   x0 = params.x0
   y0 = params.y0
   sigma_x = params.sigma_x
   sigma_y = params.sigma_y
   deltaF_x = params.deltaF_x
   deltaF_y = params.deltaF_y
   F0 = params.F0
   Theta = params.Theta
ENDELSE

chan = STRCOMPRESS(channel,/REMOVE_ALL)
ch = FIX(channel)

;;; Determine pixel X and Y phase
xphase = x - ROUND(x)   
yphase = y - ROUND(y)

NPTS = N_ELEMENTS(xphase)
IF NFLUX GT 1 THEN response = make_array(size=SIZE(observed_flux)) ELSE response = fltarr(1)
ich1 = WHERE(ch EQ 1,nch1)
ich2 = WHERE(ch EQ 2,nch2)

IF nch1 NE 0 THEN response[ich1] = pixel_phase_double_gauss(xphase[ich1],yphase[ich1],x0[0],y0[0],$
                                                                 sigma_x[0],sigma_y[0],$
                                                                 deltaF_x[0],deltaF_y[0],F0[0],theta[0]/!radeg)
IF nch2 NE 0 THEN response[ich2] = pixel_phase_double_gauss(xphase[ich2],yphase[ich2],x0[1],y0[1],$
                                                                 sigma_x[1],sigma_y[1],$
                                                                 deltaF_x[1],deltaF_y[1],F0[1],theta[1]/!radeg)

RETURN,observed_flux/response

END

