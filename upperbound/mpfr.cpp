#include <array>
#include <mpfr.h>

static const int precision = 64;
using rxystr_t = std::array<const char *, 9>;

#define CHECK(condition)    \
    do {                    \
        if (!(condition)) { \
            abort();        \
        }                   \
    } while (false)

class Mympfr
{
  public:
    static_assert((precision >= MPFR_PREC_MIN) && (precision <= MPFR_PREC_MAX), "");

    Mympfr() { mpfr_init2(value, precision); }
    ~Mympfr() { mpfr_clear(value); }

    operator mpfr_ptr() { return value; }
    mpfr_ptr operator->() { return value; }

  private:
    Mympfr(const Mympfr &) = delete;
    Mympfr &operator=(const Mympfr &) = delete;

    mpfr_t value;
};

class Verifier
{
  public:
    Verifier();

    void verify(const char *ratestr, const char *upperboundstr, const rxystr_t &rxystr);

  private:
    Mympfr rate;
    Mympfr tmpa;
    Mympfr tmpb;
    Mympfr upperbound;
    std::array<Mympfr, 3> rx;
    std::array<Mympfr, 3> ry;
    std::array<Mympfr, 9> pxy;
    std::array<Mympfr, 9> rxy;
};

Verifier::Verifier()
{
    for (uint x = 0; x < 3; ++x) {
        for (uint y = 0; y < 3; ++y) {
            CHECK(mpfr_set_si(tmpa, (x == y) ? 6 : 9997, MPFR_RNDN) == 0);
            mpfr_div_si(pxy.at(3 * x + y), tmpa, 60000, MPFR_RNDD);
        }
    }
}

void Verifier::verify(const char *ratestr, const char *upperboundstr, const rxystr_t &rxystr)
{
    // initialize values and perform basic checks

    CHECK(mpfr_strtofr(rate, ratestr, nullptr, 0, MPFR_RNDN) == 0);
    CHECK(mpfr_cmp_si(rate, 0) > 0);
    CHECK(mpfr_cmp_si(rate, 1) < 0);

    CHECK(mpfr_strtofr(upperbound, upperboundstr, nullptr, 0, MPFR_RNDN) == 0);
    CHECK(mpfr_cmp_si(upperbound, 0) > 0);
    CHECK(mpfr_cmp_si(upperbound, 1) < 0);

    for (uint i = 0; i < 9; ++i) {
        CHECK(mpfr_strtofr(rxy.at(i), rxystr.at(i), nullptr, 0, MPFR_RNDN) == 0);
        CHECK(mpfr_cmp_si(rxy.at(i), 0) > 0);
        CHECK(mpfr_cmp_si(rxy.at(i), 1) < 0);
    }

    // check that rxy is a probability mass function

    CHECK(mpfr_set_si(tmpa, -1, MPFR_RNDN) == 0);

    for (uint i = 0; i < 9; ++i) {
        CHECK(mpfr_add(tmpa, tmpa, rxy.at(i), MPFR_RNDN) == 0);
    }

    CHECK(mpfr_zero_p(tmpa) != 0);

    // compute rx and ry

    for (uint x = 0; x < 3; ++x) {
        mpfr_set_zero(rx.at(x), 0);

        for (uint y = 0; y < 3; ++y) {
            CHECK(mpfr_add(rx.at(x), rx.at(x), rxy.at(3 * x + y), MPFR_RNDN) == 0);
        }
    }

    for (uint y = 0; y < 3; ++y) {
        mpfr_set_zero(ry.at(y), 0);

        for (uint x = 0; x < 3; ++x) {
            CHECK(mpfr_add(ry.at(y), ry.at(y), rxy.at(3 * x + y), MPFR_RNDN) == 0);
        }
    }

    // check that D(rxy||rxry) < rate

    mpfr_set_zero(tmpa, 0);

    for (uint x = 0; x < 3; ++x) {
        for (uint y = 0; y < 3; ++y) {
            mpfr_mul(tmpb, rx.at(x), ry.at(y), MPFR_RNDD);
            mpfr_div(tmpb, rxy.at(3 * x + y), tmpb, MPFR_RNDU);
            mpfr_log(tmpb, tmpb, MPFR_RNDU);
            mpfr_mul(tmpb, rxy.at(3 * x + y), tmpb, MPFR_RNDU);
            mpfr_add(tmpa, tmpa, tmpb, MPFR_RNDU);
        }
    }

    CHECK(mpfr_cmp(tmpa, rate) < 0);
    mpfr_printf("%.20RUf\n", static_cast<mpfr_ptr>(tmpa));

    // check that D(rxy||pxy) < upperbound

    mpfr_set_zero(tmpa, 0);

    for (uint i = 0; i < 9; ++i) {
        mpfr_div(tmpb, rxy.at(i), pxy.at(i), MPFR_RNDU);
        mpfr_log(tmpb, tmpb, MPFR_RNDU);
        mpfr_mul(tmpb, rxy.at(i), tmpb, MPFR_RNDU);
        mpfr_add(tmpa, tmpa, tmpb, MPFR_RNDU);
    }

    CHECK(mpfr_cmp(tmpa, upperbound) < 0);
    mpfr_printf("%.20RUf\n", static_cast<mpfr_ptr>(tmpa));
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
