//
// MATLAB Compiler: 6.1 (R2015b)
// Date: Mon Sep 28 16:51:10 2015
// Arguments: "-B" "macro_default" "-W" "cpplib:covarep_cpp" "-T" "link:lib"
// "-d" "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\covarep_cpp\for_testing"
// "-v"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\creaky_voice_det
// ection\detect_creaky_voice.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\envelope\env_te.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\gci_sedreams.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\get_vq_params.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\vocoder\hmpd\hmpd_analysis.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\vocoder\hmpd\hmpd_analysis_fea
// tures.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\envelope\hspec2fwcep.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\iaif_gci.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\lpcresidual.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\mdq.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\peakslope.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\pitch_srh.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\polarity_reskew.
// m" "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\glottalsource\rd_msp.m"
// "\\bigvh\cicero\twoertwein_dcaps\covarep-1.3.2\sinusoidal\sin_analysis.m" 
//

#include <stdio.h>
#define EXPORTING_covarep_cpp 1
#include "covarep_cpp.h"

static HMCRINSTANCE _mcr_inst = NULL;


#if defined( _MSC_VER) || defined(__BORLANDC__) || defined(__WATCOMC__) || defined(__LCC__)
#ifdef __LCC__
#undef EXTERN_C
#endif
#include <windows.h>

static char path_to_dll[_MAX_PATH];

BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, void *pv)
{
    if (dwReason == DLL_PROCESS_ATTACH)
    {
        if (GetModuleFileName(hInstance, path_to_dll, _MAX_PATH) == 0)
            return FALSE;
    }
    else if (dwReason == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}
#endif
#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultPrintHandler(const char *s)
{
  return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultErrorHandler(const char *s)
{
  int written = 0;
  size_t len = 0;
  len = strlen(s);
  written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
  if (len > 0 && s[ len-1 ] != '\n')
    written += mclWrite(2 /* stderr */, "\n", sizeof(char));
  return written;
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_covarep_cpp_C_API
#define LIB_covarep_cpp_C_API /* No special import/export declaration */
#endif

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV covarep_cppInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
  if (_mcr_inst != NULL)
    return true;
  if (!mclmcrInitialize())
    return false;
  if (!GetModuleFileName(GetModuleHandle("covarep_cpp"), path_to_dll, _MAX_PATH))
    return false;
    {
        mclCtfStream ctfStream = 
            mclGetEmbeddedCtfStream(path_to_dll);
        if (ctfStream) {
            bResult = mclInitializeComponentInstanceEmbedded(   &_mcr_inst,
                                                                error_handler, 
                                                                print_handler,
                                                                ctfStream);
            mclDestroyStream(ctfStream);
        } else {
            bResult = 0;
        }
    }  
    if (!bResult)
    return false;
  return true;
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV covarep_cppInitialize(void)
{
  return covarep_cppInitializeWithHandlers(mclDefaultErrorHandler, 
                                           mclDefaultPrintHandler);
}

LIB_covarep_cpp_C_API 
void MW_CALL_CONV covarep_cppTerminate(void)
{
  if (_mcr_inst != NULL)
    mclTerminateInstance(&_mcr_inst);
}

LIB_covarep_cpp_C_API 
void MW_CALL_CONV covarep_cppPrintStackTrace(void) 
{
  char** stackTrace;
  int stackDepth = mclGetStackTrace(&stackTrace);
  int i;
  for(i=0; i<stackDepth; i++)
  {
    mclWrite(2 /* stderr */, stackTrace[i], sizeof(char)*strlen(stackTrace[i]));
    mclWrite(2 /* stderr */, "\n", sizeof(char)*strlen("\n"));
  }
  mclFreeStackTrace(&stackTrace, stackDepth);
}


LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxDetect_creaky_voice(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                         *prhs[])
{
  return mclFeval(_mcr_inst, "detect_creaky_voice", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxEnv_te(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "env_te", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxGci_sedreams(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "gci_sedreams", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxGet_vq_params(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "get_vq_params", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxHmpd_analysis(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "hmpd_analysis", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxHmpd_analysis_features(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                            *prhs[])
{
  return mclFeval(_mcr_inst, "hmpd_analysis_features", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxHspec2fwcep(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "hspec2fwcep", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxIaif_gci(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "iaif_gci", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxLpcresidual(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "lpcresidual", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxMdq(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "mdq", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxPeakslope(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "peakslope", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxPitch_srh(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "pitch_srh", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxPolarity_reskew(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "polarity_reskew", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxRd_msp(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "rd_msp", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_C_API 
bool MW_CALL_CONV mlxSin_analysis(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "sin_analysis", nlhs, plhs, nrhs, prhs);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV detect_creaky_voice(int nargout, mwArray& creak_pp, mwArray& creak_bin, 
                                      const mwArray& x, const mwArray& fs)
{
  mclcppMlfFeval(_mcr_inst, "detect_creaky_voice", nargout, 2, 2, &creak_pp, &creak_bin, &x, &fs);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV env_te(int nargout, mwArray& E, mwArray& cc, mwArray& n, const mwArray& 
                         S, const mwArray& order, const mwArray& winlen, const mwArray& 
                         mode, const mwArray& maxit, const mwArray& prec, const mwArray& 
                         presmooth_factor, const mwArray& gamma)
{
  mclcppMlfFeval(_mcr_inst, "env_te", nargout, 3, 8, &E, &cc, &n, &S, &order, &winlen, &mode, &maxit, &prec, &presmooth_factor, &gamma);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV gci_sedreams(int nargout, mwArray& gci, mwArray& MeanBasedSignal, 
                               mwArray& res, const mwArray& wave, const mwArray& fs, 
                               const mwArray& f0mean, const mwArray& polarity, const 
                               mwArray& opt)
{
  mclcppMlfFeval(_mcr_inst, "gci_sedreams", nargout, 3, 5, &gci, &MeanBasedSignal, &res, &wave, &fs, &f0mean, &polarity, &opt);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV get_vq_params(int nargout, mwArray& NAQ, mwArray& QOQ, mwArray& H1H2, 
                                mwArray& HRF, mwArray& PSP, const mwArray& gf, const 
                                mwArray& gfd, const mwArray& fs, const mwArray& GCI)
{
  mclcppMlfFeval(_mcr_inst, "get_vq_params", nargout, 5, 4, &NAQ, &QOQ, &H1H2, &HRF, &PSP, &gf, &gfd, &fs, &GCI);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV hmpd_analysis(int nargout, mwArray& f0s, mwArray& AE, mwArray& PDM, 
                                mwArray& PDD, mwArray& opt, const mwArray& wav, const 
                                mwArray& fs, const mwArray& f0s_in1, const mwArray& 
                                opt_in1)
{
  mclcppMlfFeval(_mcr_inst, "hmpd_analysis", nargout, 5, 4, &f0s, &AE, &PDM, &PDD, &opt, &wav, &fs, &f0s_in1, &opt_in1);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV hmpd_analysis_features(int nargout, mwArray& f0s, mwArray& AE, mwArray& 
                                         PDM, mwArray& PDD, mwArray& opt, const mwArray& 
                                         frames, const mwArray& fs, const mwArray& 
                                         opt_in1)
{
  mclcppMlfFeval(_mcr_inst, "hmpd_analysis_features", nargout, 5, 3, &f0s, &AE, &PDM, &PDD, &opt, &frames, &fs, &opt_in1);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV hspec2fwcep(int nargout, mwArray& fwcep, const mwArray& C, const 
                              mwArray& fs, const mwArray& order, const mwArray& warpfn, 
                              const mwArray& varargin)
{
  mclcppMlfFeval(_mcr_inst, "hspec2fwcep", nargout, 1, -5, &fwcep, &C, &fs, &order, &warpfn, &varargin);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV iaif_gci(int nargout, mwArray& g, mwArray& gd, mwArray& a, mwArray& ag, 
                           const mwArray& x, const mwArray& fs, const mwArray& GCI, const 
                           mwArray& p_vt, const mwArray& p_gl, const mwArray& d, const 
                           mwArray& hpfilt)
{
  mclcppMlfFeval(_mcr_inst, "iaif_gci", nargout, 4, 7, &g, &gd, &a, &ag, &x, &fs, &GCI, &p_vt, &p_gl, &d, &hpfilt);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV lpcresidual(int nargout, mwArray& res, mwArray& LPCcoeff, const 
                              mwArray& x, const mwArray& L, const mwArray& shift, const 
                              mwArray& order)
{
  mclcppMlfFeval(_mcr_inst, "lpcresidual", nargout, 2, 4, &res, &LPCcoeff, &x, &L, &shift, &order);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV mdq(int nargout, mwArray& m, const mwArray& res, const mwArray& fs, 
                      const mwArray& GCI)
{
  mclcppMlfFeval(_mcr_inst, "mdq", nargout, 1, 3, &m, &res, &fs, &GCI);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV peakslope(int nargout, mwArray& PS, const mwArray& s, const mwArray& fs)
{
  mclcppMlfFeval(_mcr_inst, "peakslope", nargout, 1, 2, &PS, &s, &fs);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV pitch_srh(int nargout, mwArray& F0s, mwArray& VUVDecisions, mwArray& 
                            SRHVal, mwArray& time, const mwArray& wave, const mwArray& 
                            fs, const mwArray& f0min, const mwArray& f0max, const 
                            mwArray& hopsize)
{
  mclcppMlfFeval(_mcr_inst, "pitch_srh", nargout, 4, 5, &F0s, &VUVDecisions, &SRHVal, &time, &wave, &fs, &f0min, &f0max, &hopsize);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV polarity_reskew(int nargout, mwArray& polarity, const mwArray& s, const 
                                  mwArray& fs, const mwArray& opt)
{
  mclcppMlfFeval(_mcr_inst, "polarity_reskew", nargout, 1, 3, &polarity, &s, &fs, &opt);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV rd_msp(int nargout, mwArray& rds, const mwArray& frames, const mwArray& 
                         fs, const mwArray& opt)
{
  mclcppMlfFeval(_mcr_inst, "rd_msp", nargout, 1, 3, &rds, &frames, &fs, &opt);
}

LIB_covarep_cpp_CPP_API 
void MW_CALL_CONV sin_analysis(int nargout, mwArray& frames, mwArray& syn, mwArray& opt, 
                               const mwArray& wav, const mwArray& fs, const mwArray& f0s, 
                               const mwArray& opt_in1)
{
  mclcppMlfFeval(_mcr_inst, "sin_analysis", nargout, 3, 4, &frames, &syn, &opt, &wav, &fs, &f0s, &opt_in1);
}

