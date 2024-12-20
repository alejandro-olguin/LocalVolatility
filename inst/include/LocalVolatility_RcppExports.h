// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_LocalVolatility_RCPPEXPORTS_H_GEN_
#define RCPP_LocalVolatility_RCPPEXPORTS_H_GEN_

#include <Rcpp.h>

namespace LocalVolatility {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("LocalVolatility", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("LocalVolatility", "_LocalVolatility_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in LocalVolatility");
            }
        }
    }

    inline double american_option_lv(double s_0, double k, double tau, double r_d, double q, NumericMatrix sigma, String type, double s_min, double s_max, int n_s, int n_t, double lambda, double tolerance) {
        typedef SEXP(*Ptr_american_option_lv)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_american_option_lv p_american_option_lv = NULL;
        if (p_american_option_lv == NULL) {
            validateSignature("double(*american_option_lv)(double,double,double,double,double,NumericMatrix,String,double,double,int,int,double,double)");
            p_american_option_lv = (Ptr_american_option_lv)R_GetCCallable("LocalVolatility", "_LocalVolatility_american_option_lv");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_american_option_lv(Shield<SEXP>(Rcpp::wrap(s_0)), Shield<SEXP>(Rcpp::wrap(k)), Shield<SEXP>(Rcpp::wrap(tau)), Shield<SEXP>(Rcpp::wrap(r_d)), Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(type)), Shield<SEXP>(Rcpp::wrap(s_min)), Shield<SEXP>(Rcpp::wrap(s_max)), Shield<SEXP>(Rcpp::wrap(n_s)), Shield<SEXP>(Rcpp::wrap(n_t)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(tolerance)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double european_option_lv(double s_0, double k, double tau, double r_d, double q, NumericMatrix sigma, String type, double s_min, double s_max, int n_s, int n_t) {
        typedef SEXP(*Ptr_european_option_lv)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_european_option_lv p_european_option_lv = NULL;
        if (p_european_option_lv == NULL) {
            validateSignature("double(*european_option_lv)(double,double,double,double,double,NumericMatrix,String,double,double,int,int)");
            p_european_option_lv = (Ptr_european_option_lv)R_GetCCallable("LocalVolatility", "_LocalVolatility_european_option_lv");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_european_option_lv(Shield<SEXP>(Rcpp::wrap(s_0)), Shield<SEXP>(Rcpp::wrap(k)), Shield<SEXP>(Rcpp::wrap(tau)), Shield<SEXP>(Rcpp::wrap(r_d)), Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(type)), Shield<SEXP>(Rcpp::wrap(s_min)), Shield<SEXP>(Rcpp::wrap(s_max)), Shield<SEXP>(Rcpp::wrap(n_s)), Shield<SEXP>(Rcpp::wrap(n_t)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double american_option_lv_2d(double s_0, double x_0, double k, double tau, double r_d, double r_f, double q, NumericMatrix sigma_s, NumericMatrix sigma_x, double rho, String type, double s_min, double s_max, double x_min, double x_max, int n_s, int n_x, int n_t, double alpha, double lambda, double tolerance) {
        typedef SEXP(*Ptr_american_option_lv_2d)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_american_option_lv_2d p_american_option_lv_2d = NULL;
        if (p_american_option_lv_2d == NULL) {
            validateSignature("double(*american_option_lv_2d)(double,double,double,double,double,double,double,NumericMatrix,NumericMatrix,double,String,double,double,double,double,int,int,int,double,double,double)");
            p_american_option_lv_2d = (Ptr_american_option_lv_2d)R_GetCCallable("LocalVolatility", "_LocalVolatility_american_option_lv_2d");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_american_option_lv_2d(Shield<SEXP>(Rcpp::wrap(s_0)), Shield<SEXP>(Rcpp::wrap(x_0)), Shield<SEXP>(Rcpp::wrap(k)), Shield<SEXP>(Rcpp::wrap(tau)), Shield<SEXP>(Rcpp::wrap(r_d)), Shield<SEXP>(Rcpp::wrap(r_f)), Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(sigma_s)), Shield<SEXP>(Rcpp::wrap(sigma_x)), Shield<SEXP>(Rcpp::wrap(rho)), Shield<SEXP>(Rcpp::wrap(type)), Shield<SEXP>(Rcpp::wrap(s_min)), Shield<SEXP>(Rcpp::wrap(s_max)), Shield<SEXP>(Rcpp::wrap(x_min)), Shield<SEXP>(Rcpp::wrap(x_max)), Shield<SEXP>(Rcpp::wrap(n_s)), Shield<SEXP>(Rcpp::wrap(n_x)), Shield<SEXP>(Rcpp::wrap(n_t)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(tolerance)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double american_option_2d(double s_0, double x_0, double k, double tau, double r_d, double r_f, double q, double sigma_s, double sigma_x, double rho, String type, double s_min, double s_max, double x_min, double x_max, int n_s, int n_x, int n_t, double alpha, double lambda, double tolerance) {
        typedef SEXP(*Ptr_american_option_2d)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_american_option_2d p_american_option_2d = NULL;
        if (p_american_option_2d == NULL) {
            validateSignature("double(*american_option_2d)(double,double,double,double,double,double,double,double,double,double,String,double,double,double,double,int,int,int,double,double,double)");
            p_american_option_2d = (Ptr_american_option_2d)R_GetCCallable("LocalVolatility", "_LocalVolatility_american_option_2d");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_american_option_2d(Shield<SEXP>(Rcpp::wrap(s_0)), Shield<SEXP>(Rcpp::wrap(x_0)), Shield<SEXP>(Rcpp::wrap(k)), Shield<SEXP>(Rcpp::wrap(tau)), Shield<SEXP>(Rcpp::wrap(r_d)), Shield<SEXP>(Rcpp::wrap(r_f)), Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(sigma_s)), Shield<SEXP>(Rcpp::wrap(sigma_x)), Shield<SEXP>(Rcpp::wrap(rho)), Shield<SEXP>(Rcpp::wrap(type)), Shield<SEXP>(Rcpp::wrap(s_min)), Shield<SEXP>(Rcpp::wrap(s_max)), Shield<SEXP>(Rcpp::wrap(x_min)), Shield<SEXP>(Rcpp::wrap(x_max)), Shield<SEXP>(Rcpp::wrap(n_s)), Shield<SEXP>(Rcpp::wrap(n_x)), Shield<SEXP>(Rcpp::wrap(n_t)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(tolerance)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double european_option_lv_2d(double s_0, double x_0, double k, double tau, double r_d, double r_f, double q, NumericMatrix sigma_s, NumericMatrix sigma_x, double rho, String type, double s_min, double s_max, double x_min, double x_max, int n_s, int n_x, int n_t, double alpha) {
        typedef SEXP(*Ptr_european_option_lv_2d)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_european_option_lv_2d p_european_option_lv_2d = NULL;
        if (p_european_option_lv_2d == NULL) {
            validateSignature("double(*european_option_lv_2d)(double,double,double,double,double,double,double,NumericMatrix,NumericMatrix,double,String,double,double,double,double,int,int,int,double)");
            p_european_option_lv_2d = (Ptr_european_option_lv_2d)R_GetCCallable("LocalVolatility", "_LocalVolatility_european_option_lv_2d");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_european_option_lv_2d(Shield<SEXP>(Rcpp::wrap(s_0)), Shield<SEXP>(Rcpp::wrap(x_0)), Shield<SEXP>(Rcpp::wrap(k)), Shield<SEXP>(Rcpp::wrap(tau)), Shield<SEXP>(Rcpp::wrap(r_d)), Shield<SEXP>(Rcpp::wrap(r_f)), Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(sigma_s)), Shield<SEXP>(Rcpp::wrap(sigma_x)), Shield<SEXP>(Rcpp::wrap(rho)), Shield<SEXP>(Rcpp::wrap(type)), Shield<SEXP>(Rcpp::wrap(s_min)), Shield<SEXP>(Rcpp::wrap(s_max)), Shield<SEXP>(Rcpp::wrap(x_min)), Shield<SEXP>(Rcpp::wrap(x_max)), Shield<SEXP>(Rcpp::wrap(n_s)), Shield<SEXP>(Rcpp::wrap(n_x)), Shield<SEXP>(Rcpp::wrap(n_t)), Shield<SEXP>(Rcpp::wrap(alpha)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double european_option_2d(double s0, double x0, double k, double tau, double r_d, double r_f, double q, double sigma_s, double sigma_x, double rho, String type, double s_min, double s_max, double x_min, double x_max, int n_s, int n_x, int n_t, double alpha) {
        typedef SEXP(*Ptr_european_option_2d)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_european_option_2d p_european_option_2d = NULL;
        if (p_european_option_2d == NULL) {
            validateSignature("double(*european_option_2d)(double,double,double,double,double,double,double,double,double,double,String,double,double,double,double,int,int,int,double)");
            p_european_option_2d = (Ptr_european_option_2d)R_GetCCallable("LocalVolatility", "_LocalVolatility_european_option_2d");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_european_option_2d(Shield<SEXP>(Rcpp::wrap(s0)), Shield<SEXP>(Rcpp::wrap(x0)), Shield<SEXP>(Rcpp::wrap(k)), Shield<SEXP>(Rcpp::wrap(tau)), Shield<SEXP>(Rcpp::wrap(r_d)), Shield<SEXP>(Rcpp::wrap(r_f)), Shield<SEXP>(Rcpp::wrap(q)), Shield<SEXP>(Rcpp::wrap(sigma_s)), Shield<SEXP>(Rcpp::wrap(sigma_x)), Shield<SEXP>(Rcpp::wrap(rho)), Shield<SEXP>(Rcpp::wrap(type)), Shield<SEXP>(Rcpp::wrap(s_min)), Shield<SEXP>(Rcpp::wrap(s_max)), Shield<SEXP>(Rcpp::wrap(x_min)), Shield<SEXP>(Rcpp::wrap(x_max)), Shield<SEXP>(Rcpp::wrap(n_s)), Shield<SEXP>(Rcpp::wrap(n_x)), Shield<SEXP>(Rcpp::wrap(n_t)), Shield<SEXP>(Rcpp::wrap(alpha)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

}

#endif // RCPP_LocalVolatility_RCPPEXPORTS_H_GEN_
