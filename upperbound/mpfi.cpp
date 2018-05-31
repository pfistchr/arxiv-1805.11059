#include <array>
#include <mpfi.h>

static const int precision = 64;
using rxystr_t = std::array<const char *, 9>;

#define CHECK(condition)    \
    do {                    \
        if (!(condition)) { \
            abort();        \
        }                   \
    } while (false)

class Mympfi
{
  public:
    static_assert((precision >= MPFR_PREC_MIN) && (precision <= MPFR_PREC_MAX), "");

    Mympfi() { mpfi_init2(value, precision); }
    ~Mympfi() { mpfi_clear(value); }

    operator mpfi_ptr() { return value; }
    mpfi_ptr operator->() { return value; }

  private:
    Mympfi(const Mympfi &) = delete;
    Mympfi &operator=(const Mympfi &) = delete;

    mpfi_t value;
};

class Verifier
{
  public:
    Verifier();

    void verify(const char *ratestr, const char *upperboundstr, const rxystr_t &rxystr);

  private:
    Mympfi rate;
    Mympfi tmpa;
    Mympfi tmpb;
    Mympfi upperbound;
    std::array<Mympfi, 3> rx;
    std::array<Mympfi, 3> ry;
    std::array<Mympfi, 9> pxy;
    std::array<Mympfi, 9> rxy;
};

Verifier::Verifier()
{
    for (uint x = 0; x < 3; ++x) {
        for (uint y = 0; y < 3; ++y) {
            mpfi_set_si(tmpa, (x == y) ? 6 : 9997);
            mpfi_div_si(pxy.at(3 * x + y), tmpa, 60000);
        }
    }
}

void Verifier::verify(const char *ratestr, const char *upperboundstr, const rxystr_t &rxystr)
{
    // initialize values and perform basic checks

    mpfi_set_str(rate, ratestr, 0);
    CHECK(mpfi_cmp_si(rate, 0) > 0);
    CHECK(mpfi_cmp_si(rate, 1) < 0);

    mpfi_set_str(upperbound, upperboundstr, 0);
    CHECK(mpfi_cmp_si(upperbound, 0) > 0);
    CHECK(mpfi_cmp_si(upperbound, 1) < 0);

    for (uint i = 0; i < 9; ++i) {
        mpfi_set_str(rxy.at(i), rxystr.at(i), 0);
        CHECK(mpfi_cmp_si(rxy.at(i), 0) > 0);
        CHECK(mpfi_cmp_si(rxy.at(i), 1) < 0);
    }

    // check that rxy is a probability mass function

    mpfi_set_si(tmpa, -1);

    for (uint i = 0; i < 9; ++i) {
        mpfi_add(tmpa, tmpa, rxy.at(i));
    }

    CHECK(mpfi_is_zero(tmpa) != 0);

    // compute rx and ry

    for (uint x = 0; x < 3; ++x) {
        mpfi_set_si(rx.at(x), 0);

        for (uint y = 0; y < 3; ++y) {
            mpfi_add(rx.at(x), rx.at(x), rxy.at(3 * x + y));
        }
    }

    for (uint y = 0; y < 3; ++y) {
        mpfi_set_si(ry.at(y), 0);

        for (uint x = 0; x < 3; ++x) {
            mpfi_add(ry.at(y), ry.at(y), rxy.at(3 * x + y));
        }
    }

    // check that D(rxy||rxry) < rate

    mpfi_set_si(tmpa, 0);

    for (uint x = 0; x < 3; ++x) {
        for (uint y = 0; y < 3; ++y) {
            mpfi_mul(tmpb, rx.at(x), ry.at(y));
            mpfi_div(tmpb, rxy.at(3 * x + y), tmpb);
            mpfi_log(tmpb, tmpb);
            mpfi_mul(tmpb, rxy.at(3 * x + y), tmpb);
            mpfi_add(tmpa, tmpa, tmpb);
        }
    }

    CHECK(mpfi_cmp(tmpa, rate) < 0);
    mpfr_printf("%.20RUf\n", static_cast<mpfr_ptr>(&tmpa->right));

    // check that D(rxy||pxy) < upperbound

    mpfi_set_si(tmpa, 0);

    for (uint i = 0; i < 9; ++i) {
        mpfi_div(tmpb, rxy.at(i), pxy.at(i));
        mpfi_log(tmpb, tmpb);
        mpfi_mul(tmpb, rxy.at(i), tmpb);
        mpfi_add(tmpa, tmpa, tmpb);
    }

    CHECK(mpfi_cmp(tmpa, upperbound) < 0);
    mpfr_printf("%.20RUf\n", static_cast<mpfr_ptr>(&tmpa->right));
}

static const char *rateA = "0x0.079d"; // 3898 / 2^17
static const char *rateB = "0x0.07c8"; // 3984 / 2^17

static const char *upperboundA = "0x0.d02a7208f52317"; // 58593464420737815 / 2^56
static const char *upperboundB = "0x0.cf6aa03d0f0253"; // 58382556630811219 / 2^56

static const rxystr_t rxyA = {"0x0.00d18e2d53dba4", "0x0.6c6ebcb6c6ea40", "0x0.6c6ebcb6c6ea40",
                              "0x0.006ff71d804e2a", "0x0.03d405476786bd", "0x0.0ee47fcda75307",
                              "0x0.006ff71d804e2a", "0x0.0ee47fcda75307", "0x0.03d405476786bd"};

static const rxystr_t rxyB = {"0x0.0184ae0a6be14a", "0x0.35ba25f4e1fd7f", "0x0.870eb8aa072ec5",
                              "0x0.02598735ff8940", "0x0.057e6f74c876f3", "0x0.35ba25f4e1fd7f",
                              "0x0.00422176958a36", "0x0.02598735ff8940", "0x0.0184ae0a6be14a"};

int main()
{
    Verifier verifier;
    verifier.verify(rateA, upperboundA, rxyA);
    verifier.verify(rateB, upperboundB, rxyB);
    printf("finish\n");
    return 0;
}
