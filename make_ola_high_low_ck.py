"""
Python script using module spiceypy to convert time-invariant
fixed-offset (class 4) reference frame (reffrm) in FK to time-dependent
CK-based (class 3) reference frame in [FK + CK]

Usage:

  [DEBUG=] [VERBOSE=] python make_ola_high_low_ck.py [--create] [--test]

Associated sample input files:

  make_ola_high_low_ck.py  - this script; also a meta-kernel
  naif0012.tls             - LSK, valid ca. late 2018
  ORX_SCLKSCET.00039.tsc   - SCLK, OSIRIS-REx (ORX), ca. late 2018
  orx_v11.tf               - ORX FK version 1.1, pre-2018-10-25
  orx_v12.tf               - ORX FK version 1.2, post-2018-10.25
  orx_v13_draft.tf         - ORX FK if selected frames become class 3

Associated sample output file:

  orx_ola_fk_replacement.bc  - CK to work with orx_v13_draft.tf

"""
import os
import sys
import pprint
import spiceypy as sp
import traceback as tb

"""
Sample meta-kernel (MK) to control this script; contains five Kernel
Pool Variables (KPVs):

1) KERNELS_TO_LOAD - for FURNSH/KEEPER susbsystem; LSK, SCLK, new FK)
2) FRAMES_TO_CONVERT - names of frames that change at various epochs
3) FKS - names of old FKs that are valid over known windows
4) STOP_ETS - the window endpoints; the last is -1d32
5) OUTPUT_CK - the CK to be written to replace the old FKs' transforms

\begindata
KERNELS_TO_LOAD = (
  'naif0012.tls'
  'orx_v13_draft.tf'
  'ORX_SCLKSCET.00039.tsc'
)

FRAMES_TO_CONVERT = (
  'ORX_OLA_HIGH'
  'ORX_OLA_LOW'
)

FKS = (
  'orx_v11.tf'
  'orx_v12.tf'
)

STOP_ETS = (
  @2018-10-25/00:01:09.184
  -1d32
)

OUTPUT_CK = (
  'orx_ola_fk_replacement.bc'
)

\begintext
"""

### Debugging and logging
doDebug = "DEBUG" in os.environ
doVerbose = "VERBOSE" in os.environ


########################################################################
def create_ck(kernels,frames_to_convert,fks,stop_ets,output_ck):
  """
  Create CK from FKs and times at which the FKs end being valid

  """

  ### FURNSH any kernels
  for kernel in kernels: sp.furnsh(kernel)

  ### Set last ET to None so it will be initialized in loop's first pass
  last_et = None

  ### Do not overwrite existing output CK
  assert (not sp.exists(output_ck)) or dict()['Cannot overwrite existing CK[{}]'.format(output_ck)]

  ### Initialize CK handle to None
  ck_handle = None

  ### Loop over CKs
  while fks:

    ### Get FK and corresponding stop ET
    fk = fks.pop()
    stop_et = stop_ets.pop()

    ### Loop over reference frames (reffrms)
    for reffrm in frames_to_convert:

      ### Get reffrm ID, SCLK ID; N.B. latter comes from new FK
      reffrm_id = sp.gipool('FRAME_{}'.format(reffrm.upper()),0,1)[0]
      sclk_id = sp.gipool('CK_{}_{}'.format(reffrm_id,'SCLK'),0,1)[0]

      ### Set start DP-SCLK of window for this old FK (outer loop)
      ### - Set to zero for first window
      ### - Covnert ET to DP-SCLK for subsequent windows
      if last_et is None: begtim = 0.0
      else              : begtim = sp.sce2t(sclk_id,last_et)

      ### Load old FK, get RELATIVE frame name, get time-invariant
      ### matrix, and unload old FK
      sp.furnsh(fk)
      relative_reffrm = sp.gcpool('TKFRAME_{}_{}'.format(reffrm_id,'RELATIVE'),0,1,99)[0]
      mtx = sp.pxform(relative_reffrm,reffrm,0.0)
      sp.unload(fk)

      ### Covnert matrix to quaternion
      quat = sp.m2q(mtx)

      ### Calculate tick rate:  seconds per tick
      rate = (sp.sct2e(sclk_id,1e3) - sp.sct2e(sclk_id,0.)) / 1e3

      if doVerbose:
        ### if VERBOSE environment variable is present, log information
        print((relative_reffrm,reffrm,fk,last_et,stop_et,'{:010.3f}'.format(1./rate),quat,))

      ### Set stop DP-SCLK of window
      if stop_et < -1e30:
        ### Use end of encoded DP_SCLK for final window
        endtim = sp.gdpool('SCLK_PARTITION_END_{}'.format(-sclk_id),0,999)[-1]
      else:
        ### Else convert stop ET to DP-SCLK
        endtim = sp.sce2t(sclk_id,stop_et)

      if doDebug:
        ### Debug output
        pprint.pprint(dict(fk=fk
                          ,reffrm=reffrm
                          ,reffrm_id=reffrm_id
                          ,relative_reffrm=relative_reffrm
                          ,rate=rate
                          ,begtim=begtim
                          ,endtim=endtim
                          ,diff=endtim-begtim
                          ,mtx=mtx
                          ,quat=quat
                          ))


      ### Open CK, once
      if ck_handle is None: ck_handle = sp.ckopn(output_ck,'ORX FK REPLACEMENT',0)

      ### Write Type 2 CK segment with one record; angular velocity = 0
      sp.ckw02(ck_handle,begtim,endtim
              ,reffrm_id,relative_reffrm
              ,'{}[{}]'.format(os.path.basename(fk),reffrm)[:40]
              ,1,[begtim],[endtim]
              ,[quat],[[0.,0.,0.]],[rate]
              )

    ### Save stop ET for start ET of next pass
    last_et = stop_et

  ### Close CK
  if not (ck_handle is None): sp.ckcls(ck_handle)


########################################################################
def test_ck(kernels,frames_to_convert,fks,stop_ets,output_ck):
  """
  Text CK against FKs at times at which the FKs are valid

  """

  ### Load base kernels (LSK, new FK, SCLK)
  for kernel in kernels+[output_ck]: sp.furnsh(kernel)

  ### Set last ET to None so it will be initialized in loop's first pass
  last_et = None

  ### Create dict of info; keys will be new FK filenames
  dt = dict()

  while fks:

    ### Pop FK and stop ET off of lists
    fk = fks.pop()
    stop_et = stop_ets.pop()

    ### Create dict for this FK
    dt[fk] = dict()

    ### Loop over refernce frames
    for reffrm in frames_to_convert:

      ### Get reffrm ID, SCLK ID; N.B. latter comes from new FK
      reffrm_id = sp.gipool('FRAME_{}'.format(reffrm.upper()),0,1)[0]
      sclk_id = sp.gipool('CK_{}_{}'.format(reffrm_id,'SCLK'),0,1)[0]

      ### Set start DP-SCLK of window for this old FK (outer loop)
      ### - Set to zero for first window
      ### - Covnert ET to DP-SCLK for subsequent windows
      if last_et is None: et_lo = sp.sct2e(sclk_id,0.)
      else              : et_lo = last_et

      ### Load old FK, get RELATIVE frame name, get time-invariant
      ### matrix, and unload old FK
      sp.furnsh(fk)
      relative_reffrm = sp.gcpool('TKFRAME_{}_{}'.format(reffrm_id,'RELATIVE'),0,1,99)[0]
      sp.unload(fk)

      ### Get ETs at which to do the tests:
      ### - 10s after start of window
      ### - 10s before end of window, or 1e6s after start if last window
      if stop_et < -1e30:
        et_test_lo = et_lo + 10.
        et_test_hi = et_lo + 1e6
      else:
        et_delta = min([10.,(stop_et-et_lo)/3.])
        et_test_lo = et_lo + et_delta
        et_test_hi = stop_et - et_delta

      ### Save the relative reffrm, the reffrm, the window, and an empty
      ### dict for this reffrm under this FK
      dt[fk][reffrm] = (relative_reffrm,et_test_lo,et_test_hi,dict())

    ### For next pass
    last_et = stop_et

  ### Clear all kernels, and test
  sp.kclear()
  assert 0 == sp.ktotal('all')

  ### Load base kernels including new FK and new CK
  for kernel in kernels+[output_ck]: sp.furnsh(kernel)

  ### Loop over old FKs, reffrms, and ETs
  for fk in dt:

    for reffrm in dt[fk]:

      ### Retrieve relative reffrm, ETs and quat dict
      relative_reffrm,et_test_lo,et_test_hi,dtquat = dt[fk][reffrm]

      for et in (et_test_lo,et_test_hi,):

        ### Lookup CK-based matrix, convrt to and save quat at each ET
        dtquat[et] = sp.m2q(sp.pxform(relative_reffrm,reffrm,et))

  ### Loop over the old FKs again
  for fk in dt:

    ### Clear all kernels, and test
    sp.kclear()
    assert 0 == sp.ktotal('all')

    ### Load only the old FK
    sp.furnsh(fk)

    ### Loop over reffrms, and ETs
    for reffrm in dt[fk]:

      relative_reffrm,et_test_lo,et_test_hi,dtquat = dt[fk][reffrm]

      for et in (et_test_lo,et_test_hi,):

        ### Calculate norm of difference of CK-based and FK-based quats
        quat_error = round(sp.vnorm(dtquat[et]-sp.m2q(sp.pxform(relative_reffrm,reffrm,et))),16)

        ### Output that norm as an error for each case, which norm
        ### should be zero
        print(dict(fk=fk
                  ,quat_error=quat_error
                  ,relative_reffrm=relative_reffrm
                  ,reffrm=reffrm
                  ,et='{:015.4f}'.format(et)
                  ))


########################################################################
def main(argv_arg=sys.argv[1:]):
  """
  Process oommand-line arguments into input to create_ck and test_ck 
  methods

  """

  ### Process arguments

  ### Create set of non-kernel options
  non_kernel_opts = set('--test --create'.split())

  ### All arguments which are not in that set are SPICE kernels
  kernels = [arg for arg in argv_arg if not (arg in non_kernel_opts)]

  ### - Use this script as the lone SPICE kernel if none provided
  if not kernels: kernels =[__file__[-4:-1]=='.py' and __file__[:-1] or __file__]

  ### - FURNSH those kernels
  for kernel in kernels: sp.furnsh(kernel)

  ### Get frames to convert from FK to CK form
  frames_to_convert = sp.gcpool('FRAMES_TO_CONVERT',0,999,99)
  assert frames_to_convert

  ### Get FK paths; use STPOOL in case paths are long
  fks = list()
  item = 0
  while True:
    ### - Use '+' for continuation
    ### - Read paths until exception; assume that indicates the end
    try   : fks.append(sp.stpool('FKS',item,'+',2048)[0])
    except:
       if doDebug:
         tb.print_exc()
         sys.stderr.write('\nN.B. Expected error; continuing ...\n\n')
       break
    item += 1

  ### Throw exception if no FKs were provided
  assert fks

  ### Get stop times for each FK, cast to list so .pop() is available
  stop_ets = list(sp.gdpool('STOP_ETS',0,999))

  ### Ensure there is a stop time for each FK
  assert len(fks)==len(stop_ets)

  ### Reverse order of FKs so .pop() will get items in order
  fks.reverse()
  stop_ets.reverse()

  ### Get output CK path
  output_ck = sp.stpool('OUTPUT_CK',0,'+',2048)[0]

  ### Unload all kernels
  sp.kclear()

  if '--create' in argv_arg:
    ### Call method to create CK; pass copies of fks and stop_ets
    create_ck(kernels,frames_to_convert,fks[:],stop_ets[:],output_ck)

  if '--test' in argv_arg:
    ### Call test method; pass copies of fks and stop_ets
    test_ck(kernels,frames_to_convert,fks[:],stop_ets[:],output_ck)


########################################################################
if "__main__" == __name__:
  main()
