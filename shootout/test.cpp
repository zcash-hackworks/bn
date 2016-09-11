#include "common/default_types/r1cs_ppzksnark_pp.hpp"
#include "common/profiling.hpp"

// CXXFLAGS="-fPIC -DBINARY_OUTPUT -DNO_PT_COMPRESSION=1" make lib CURVE=ALT_BN128 MULTICORE=1 NO_PROCPS=1 NO_GTEST=1 NO_DOCS=1 STATIC=1 NO_SUPERCOP=1 FEATUREFLAGS=-DMONTGOMERY_OUTPUT
// g++ -std=c++11 -O3 test.cpp -o test -Isrc/ -DMULTICORE -fopenmp -DBINARY_OUTPUT -DCURVE_ALT_BN128 -DSTATIC -L. -lsnark -lgmp -lsodium

typedef libsnark::default_r1cs_ppzksnark_pp curve_pp;
typedef libsnark::default_r1cs_ppzksnark_pp::G1_type curve_G1;
typedef libsnark::default_r1cs_ppzksnark_pp::G2_type curve_G2;
typedef libsnark::default_r1cs_ppzksnark_pp::GT_type curve_GT;
typedef libsnark::default_r1cs_ppzksnark_pp::Fp_type curve_Fr;

int main() {
    curve_pp::init_public_params();
    libsnark::inhibit_profiling_info = true;
    libsnark::inhibit_profiling_counters = true;

    curve_G1 a = curve_G1::one();
    curve_G2 b = curve_G2::one();

    curve_Fr c = curve_Fr("1901").inverse();
    curve_Fr d = curve_Fr("2344").inverse();

    curve_GT acc1 = curve_GT::one();

    for (size_t i = 0; i < 10000; i++) {
        acc1 = acc1 * curve_pp::reduced_pairing(a, b);
        a = c * a;
        b = d * b;
    }

    a = curve_G1::one();
    b = curve_G2::one();

    curve_GT acc2 = curve_GT::one();

    for (size_t i = 0; i < 10000; i++) {
        acc2 = acc2 * curve_pp::reduced_pairing(a, b);
        a = d * a;
        b = c * b;
    }

    assert(acc1 == acc2);
}
