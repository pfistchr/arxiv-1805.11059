#include <array>
#include <fstream>
#include <mpfr.h>
#include <set>
#include <stack>

static const int fixedpointexponent = 29;
static const int precision = 68;
static const int32_t fixedpointone = (1 << fixedpointexponent);
static const std::array<std::array<uint8_t, 3>, 6> permutations = {{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}};
using line_t = std::array<char, 171>;

struct qxybounds_t {
    // lower.at(0) <= Q_X(1) * 2^fixedpointexponent <= upper.at(0)
    // lower.at(1) <= Q_X(2) * 2^fixedpointexponent <= upper.at(1)
    // lower.at(2) <= Q_X(3) * 2^fixedpointexponent <= upper.at(2)
    // lower.at(3) <= Q_Y(1) * 2^fixedpointexponent <= upper.at(3)
    // lower.at(4) <= Q_Y(2) * 2^fixedpointexponent <= upper.at(4)
    // lower.at(5) <= Q_Y(3) * 2^fixedpointexponent <= upper.at(5)

    std::array<int32_t, 6> lower;
    std::array<int32_t, 6> upper;
};

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
    Verifier(const char *ratestr, const char *lowerboundstr);

    void verify(const qxybounds_t &qxybounds, const line_t &line);

  private:
    Mympfr alpha;
    Mympfr lowerbound;
    Mympfr oneminusalpha;
    Mympfr qxybetamin;
    Mympfr rate;
    Mympfr tmpa;
    Mympfr tmpb;
    std::array<Mympfr, 9> beta;
    std::array<Mympfr, 9> logpxy;
};

Verifier::Verifier(const char *ratestr, const char *lowerboundstr)
{
    CHECK(mpfr_strtofr(rate, ratestr, nullptr, 0, MPFR_RNDN) == 0);
    CHECK(mpfr_cmp_si(rate, 0) > 0);
    CHECK(mpfr_cmp_si(rate, 1) < 0);

    CHECK(mpfr_strtofr(lowerbound, lowerboundstr, nullptr, 0, MPFR_RNDN) == 0);
    CHECK(mpfr_cmp_si(lowerbound, 0) > 0);
    CHECK(mpfr_cmp_si(lowerbound, 1) < 0);

    for (uint x = 0; x < 3; ++x) {
        for (uint y = 0; y < 3; ++y) {
            CHECK(mpfr_set_si(tmpa, (x == y) ? 6 : 9997, MPFR_RNDN) == 0);
            mpfr_div_si(tmpa, tmpa, 60000, MPFR_RNDU);
            mpfr_log(logpxy.at(3 * x + y), tmpa, MPFR_RNDU);
        }
    }
}

void Verifier::verify(const qxybounds_t &qxybounds, const line_t &line)
{
    // initialize values and perform basic checks

    CHECK(mpfr_strtofr(alpha, &line.at(1), nullptr, 16, MPFR_RNDN) == 0);
    CHECK(mpfr_cmp_d(alpha, 0.001) > 0);
    CHECK(mpfr_cmp_d(alpha, 0.999) < 0);

    for (uint i = 0; i < 9; ++i) {
        CHECK(mpfr_strtofr(beta.at(i), &line.at(17 * i + 18), nullptr, 16, MPFR_RNDN) == 0);
        CHECK(mpfr_cmp_si(beta.at(i), 0) >= 0);
        CHECK(mpfr_cmp_si(beta.at(i), 9) < 0);
    }

    CHECK(mpfr_si_sub(oneminusalpha, 1, alpha, MPFR_RNDN) == 0);

    // determine corner points

    std::set<std::array<int32_t, 3>> qxcornerpoints;

    for (const std::array<uint8_t, 3> &permutation : permutations) {
        const int32_t qxlower = qxybounds.lower.at(permutation.at(0));
        const int32_t qxupper = qxybounds.upper.at(permutation.at(1));
        const int32_t qxother = (fixedpointone - qxlower - qxupper);

        CHECK(qxother >= qxybounds.lower.at(permutation.at(2)));
        CHECK(qxother <= qxybounds.upper.at(permutation.at(2)));

        std::array<int32_t, 3> qx = {};
        qx.at(permutation.at(0)) = qxlower;
        qx.at(permutation.at(1)) = qxupper;
        qx.at(permutation.at(2)) = qxother;
        qxcornerpoints.insert(qx);
    }

    std::set<std::array<int32_t, 3>> qycornerpoints;

    for (const std::array<uint8_t, 3> &permutation : permutations) {
        const int32_t qylower = qxybounds.lower.at(3 + permutation.at(0));
        const int32_t qyupper = qxybounds.upper.at(3 + permutation.at(1));
        const int32_t qyother = (fixedpointone - qylower - qyupper);

        CHECK(qyother >= qxybounds.lower.at(3 + permutation.at(2)));
        CHECK(qyother <= qxybounds.upper.at(3 + permutation.at(2)));

        std::array<int32_t, 3> qy = {};
        qy.at(permutation.at(0)) = qylower;
        qy.at(permutation.at(1)) = qyupper;
        qy.at(permutation.at(2)) = qyother;
        qycornerpoints.insert(qy);
    }

    // compute D = \min_j \sum_{x,y} Q_j(x,y)^{1-\alpha} \beta(x,y)

    mpfr_set_inf(qxybetamin, 0);

    for (const std::array<int32_t, 3> &qx : qxcornerpoints) {
        for (const std::array<int32_t, 3> &qy : qycornerpoints) {
            mpfr_set_zero(tmpa, 0);

            for (uint x = 0; x < 3; ++x) {
                for (uint y = 0; y < 3; ++y) {
                    if ((qx.at(x) == 0) || (qy.at(y) == 0)) {
                        continue;
                    }

                    CHECK(mpfr_set_si(tmpb, qx.at(x), MPFR_RNDN) == 0);
                    CHECK(mpfr_mul_si(tmpb, tmpb, qy.at(y), MPFR_RNDN) == 0);
                    CHECK(mpfr_div_2si(tmpb, tmpb, 2 * fixedpointexponent, MPFR_RNDN) == 0);
                    mpfr_log(tmpb, tmpb, MPFR_RNDD);
                    mpfr_mul(tmpb, tmpb, oneminusalpha, MPFR_RNDD);
                    CHECK(mpfr_number_p(tmpb) != 0);
                    mpfr_exp(tmpb, tmpb, MPFR_RNDD);
                    mpfr_mul(tmpb, tmpb, beta.at(3 * x + y), MPFR_RNDD);
                    mpfr_add(tmpa, tmpa, tmpb, MPFR_RNDD);
                }
            }

            CHECK(mpfr_number_p(tmpa) != 0);
            CHECK(mpfr_min(qxybetamin, qxybetamin, tmpa, MPFR_RNDN) == 0);
        }
    }

    // compute \left[ \sum_{x,y} (P(x,y)^\alpha + \beta(x,y))^\frac{1}{\alpha} \right]^\alpha

    mpfr_set_zero(tmpa, 0);

    for (uint i = 0; i < 9; ++i) {
        mpfr_mul(tmpb, logpxy.at(i), alpha, MPFR_RNDU);
        mpfr_exp(tmpb, tmpb, MPFR_RNDU);
        mpfr_add(tmpb, tmpb, beta.at(i), MPFR_RNDU);
        mpfr_log(tmpb, tmpb, MPFR_RNDU);
        mpfr_div(tmpb, tmpb, alpha, MPFR_RNDU);
        mpfr_exp(tmpb, tmpb, MPFR_RNDU);
        mpfr_add(tmpa, tmpa, tmpb, MPFR_RNDU);
    }

    mpfr_log(tmpa, tmpa, MPFR_RNDU);
    mpfr_mul(tmpa, tmpa, alpha, MPFR_RNDU);
    mpfr_exp(tmpa, tmpa, MPFR_RNDU);

    // compute value = -\frac{[\log (... - D) + (1 - \alpha) \cdot rate]}{\alpha}

    mpfr_sub(tmpa, tmpa, qxybetamin, MPFR_RNDU);
    mpfr_log(tmpa, tmpa, MPFR_RNDU);
    mpfr_mul(tmpb, oneminusalpha, rate, MPFR_RNDU);
    mpfr_add(tmpa, tmpa, tmpb, MPFR_RNDU);
    mpfr_div(tmpa, tmpa, alpha, MPFR_RNDU);
    CHECK(mpfr_neg(tmpa, tmpa, MPFR_RNDN) == 0);

    // check that value > lowerbound

    CHECK(mpfr_cmp(lowerbound, tmpa) < 0);
    mpfr_printf("%.20RDf\n", static_cast<mpfr_ptr>(tmpa));
}

int main()
{
    Verifier verifier("0x0.07b28", "0x0.cfca8923023b33"); // 3941 / 2^17 and 58488010525784883 / 2^56
    std::stack<qxybounds_t> stack;
    stack.push(qxybounds_t{{0, 0, 0, 0, 0, 0}, {fixedpointone, fixedpointone, fixedpointone, fixedpointone, fixedpointone, fixedpointone}});
    std::ifstream infile("input.txt");
    line_t line = {};

    while (!stack.empty()) {
        // get top case from stack

        qxybounds_t qxybounds = stack.top();
        stack.pop();

        // tighten upper and lower bounds

        for (uint j = 0; j < 6; ++j) {
            CHECK(qxybounds.lower.at(j) >= 0);
            CHECK(qxybounds.lower.at(j) < qxybounds.upper.at(j));
            CHECK(qxybounds.upper.at(j) <= fixedpointone);
        }

        qxybounds.lower.at(0) = std::max(qxybounds.lower.at(0), fixedpointone - qxybounds.upper.at(1) - qxybounds.upper.at(2));
        qxybounds.lower.at(1) = std::max(qxybounds.lower.at(1), fixedpointone - qxybounds.upper.at(0) - qxybounds.upper.at(2));
        qxybounds.lower.at(2) = std::max(qxybounds.lower.at(2), fixedpointone - qxybounds.upper.at(0) - qxybounds.upper.at(1));
        qxybounds.lower.at(3) = std::max(qxybounds.lower.at(3), fixedpointone - qxybounds.upper.at(4) - qxybounds.upper.at(5));
        qxybounds.lower.at(4) = std::max(qxybounds.lower.at(4), fixedpointone - qxybounds.upper.at(3) - qxybounds.upper.at(5));
        qxybounds.lower.at(5) = std::max(qxybounds.lower.at(5), fixedpointone - qxybounds.upper.at(3) - qxybounds.upper.at(4));
        qxybounds.upper.at(0) = std::min(qxybounds.upper.at(0), fixedpointone - qxybounds.lower.at(1) - qxybounds.lower.at(2));
        qxybounds.upper.at(1) = std::min(qxybounds.upper.at(1), fixedpointone - qxybounds.lower.at(0) - qxybounds.lower.at(2));
        qxybounds.upper.at(2) = std::min(qxybounds.upper.at(2), fixedpointone - qxybounds.lower.at(0) - qxybounds.lower.at(1));
        qxybounds.upper.at(3) = std::min(qxybounds.upper.at(3), fixedpointone - qxybounds.lower.at(4) - qxybounds.lower.at(5));
        qxybounds.upper.at(4) = std::min(qxybounds.upper.at(4), fixedpointone - qxybounds.lower.at(3) - qxybounds.lower.at(5));
        qxybounds.upper.at(5) = std::min(qxybounds.upper.at(5), fixedpointone - qxybounds.lower.at(3) - qxybounds.lower.at(4));

        for (uint j = 0; j < 6; ++j) {
            CHECK(qxybounds.lower.at(j) >= 0);
            CHECK(qxybounds.lower.at(j) < qxybounds.upper.at(j));
            CHECK(qxybounds.upper.at(j) <= fixedpointone);
        }

        // process next line

        infile.getline(line.begin(), line.size());

        if (line.at(0) == 'v') {
            // verify lower bound

            for (uint k = 1; k <= 10; ++k) {
                line.at(17 * k) = 0;
            }

            verifier.verify(qxybounds, line);
            continue;
        }

        // split case

        CHECK((line.at(0) >= 'a') && (line.at(0) <= 'f'));
        const uint splitindex = uint(line.at(0) - 'a');

        // prepare split

        CHECK((qxybounds.lower.at(splitindex) % 2) == 0);
        CHECK((qxybounds.upper.at(splitindex) % 2) == 0);
        const int32_t middlevalue = ((qxybounds.lower.at(splitindex) / 2) + (qxybounds.upper.at(splitindex) / 2));

        // construct lower and upper part

        qxybounds_t lower = qxybounds;
        qxybounds_t upper = qxybounds;
        lower.upper.at(splitindex) = middlevalue;
        upper.lower.at(splitindex) = middlevalue;

        // push parts to stack

        stack.push(upper);
        stack.push(lower);
    }

    printf("finish\n");
    return 0;
}
