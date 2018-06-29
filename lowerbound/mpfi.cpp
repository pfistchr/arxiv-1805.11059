#include <array>
#include <fstream>
#include <mpfi.h>
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
    Verifier(const char *ratestr, const char *lowerboundstr);

    void verify(const qxybounds_t &qxybounds, const line_t &line);

  private:
    Mympfi alpha;
    Mympfi lowerbound;
    Mympfi oneminusalpha;
    Mympfi qxybetamin;
    Mympfi rate;
    Mympfi tmpa;
    Mympfi tmpb;
    std::array<Mympfi, 9> beta;
    std::array<Mympfi, 9> logpxy;
};

Verifier::Verifier(const char *ratestr, const char *lowerboundstr)
{
    mpfi_set_str(rate, ratestr, 0);
    CHECK(mpfi_cmp_si(rate, 0) > 0);
    CHECK(mpfi_cmp_si(rate, 1) < 0);

    mpfi_set_str(lowerbound, lowerboundstr, 0);
    CHECK(mpfi_cmp_si(lowerbound, 0) > 0);
    CHECK(mpfi_cmp_si(lowerbound, 1) < 0);

    for (uint x = 0; x < 3; ++x) {
        for (uint y = 0; y < 3; ++y) {
            mpfi_set_si(tmpa, (x == y) ? 6 : 9997);
            mpfi_div_si(tmpa, tmpa, 60000);
            mpfi_log(logpxy.at(3 * x + y), tmpa);
        }
    }
}

void Verifier::verify(const qxybounds_t &qxybounds, const line_t &line)
{
    // initialize values and perform basic checks

    mpfi_set_str(alpha, &line.at(1), 16);
    CHECK(mpfi_cmp_d(alpha, 0.001) > 0);
    CHECK(mpfi_cmp_d(alpha, 0.999) < 0);

    for (uint i = 0; i < 9; ++i) {
        mpfi_set_str(beta.at(i), &line.at(17 * i + 18), 16);
        CHECK(mpfi_cmp_si(beta.at(i), 0) >= 0);
        CHECK(mpfi_cmp_si(beta.at(i), 9) < 0);
    }

    mpfi_si_sub(oneminusalpha, 1, alpha);

    // determine extreme points

    std::set<std::array<int32_t, 3>> qxextremepoints;

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
        qxextremepoints.insert(qx);
    }

    std::set<std::array<int32_t, 3>> qyextremepoints;

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
        qyextremepoints.insert(qy);
    }

    // compute D = qxybetamin = \min_j \sum_{x,y} Q_j(x,y)^{1-\alpha} \beta(x,y)

    // [+\infty,+\infty] is not a valid interval
    mpfr_set_inf(&qxybetamin->left, 0);
    mpfr_set_inf(&qxybetamin->right, 0);

    for (const std::array<int32_t, 3> &qx : qxextremepoints) {
        for (const std::array<int32_t, 3> &qy : qyextremepoints) {
            mpfi_set_si(tmpa, 0);

            for (uint x = 0; x < 3; ++x) {
                for (uint y = 0; y < 3; ++y) {
                    if ((qx.at(x) == 0) || (qy.at(y) == 0)) {
                        continue;
                    }

                    mpfi_set_si(tmpb, qx.at(x));
                    mpfi_mul_si(tmpb, tmpb, qy.at(y));
                    mpfi_div_2si(tmpb, tmpb, 2 * fixedpointexponent);
                    mpfi_log(tmpb, tmpb);
                    mpfi_mul(tmpb, tmpb, oneminusalpha);
                    CHECK(mpfi_bounded_p(tmpb) != 0);
                    mpfi_exp(tmpb, tmpb);
                    mpfi_mul(tmpb, tmpb, beta.at(3 * x + y));
                    mpfi_add(tmpa, tmpa, tmpb);
                }
            }

            CHECK(mpfi_bounded_p(tmpa) != 0);
            CHECK(mpfr_min(&qxybetamin->left, &qxybetamin->left, &tmpa->left, MPFR_RNDN) == 0);
            CHECK(mpfr_min(&qxybetamin->right, &qxybetamin->right, &tmpa->right, MPFR_RNDN) == 0);
        }
    }

    // check that qxybetamin is a valid interval
    CHECK(mpfr_regular_p(&qxybetamin->left) != 0);
    CHECK(mpfr_regular_p(&qxybetamin->right) != 0);
    CHECK(mpfr_cmp(&qxybetamin->left, &qxybetamin->right) < 0);

    // compute \left[ \sum_{x,y} (P(x,y)^\alpha + \beta(x,y))^\frac{1}{\alpha} \right]^\alpha

    mpfi_set_si(tmpa, 0);

    for (uint i = 0; i < 9; ++i) {
        mpfi_mul(tmpb, logpxy.at(i), alpha);
        mpfi_exp(tmpb, tmpb);
        mpfi_add(tmpb, tmpb, beta.at(i));
        mpfi_log(tmpb, tmpb);
        mpfi_div(tmpb, tmpb, alpha);
        mpfi_exp(tmpb, tmpb);
        mpfi_add(tmpa, tmpa, tmpb);
    }

    mpfi_log(tmpa, tmpa);
    mpfi_mul(tmpa, tmpa, alpha);
    mpfi_exp(tmpa, tmpa);

    // compute value = -\frac{\log \{[...]^\alpha - D\} + (1 - \alpha) \cdot rate}{\alpha}

    mpfi_sub(tmpa, tmpa, qxybetamin);
    mpfi_log(tmpa, tmpa);
    mpfi_mul(tmpb, oneminusalpha, rate);
    mpfi_add(tmpa, tmpa, tmpb);
    mpfi_div(tmpa, tmpa, alpha);
    mpfi_neg(tmpa, tmpa);

    // check that value > lowerbound

    CHECK(mpfi_cmp(lowerbound, tmpa) < 0);
    mpfr_printf("%.20RDf\n", static_cast<mpfr_ptr>(&tmpa->left));
}

int main()
{
    Verifier verifier("0x0.07b28", "0x0.cfca8923023b33"); // 3941 / 2^17 and 58488010525784883 / 2^56
    std::stack<qxybounds_t> stack;
    stack.push(qxybounds_t{{0, 0, 0, 0, 0, 0}, {fixedpointone, fixedpointone, fixedpointone, fixedpointone, fixedpointone, fixedpointone}});
    std::ifstream infile("input.txt");
    line_t line = {};

    if (!infile) {
        printf("cannot open input.txt\n");
        return 1;
    }

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
