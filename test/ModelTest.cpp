#include "gtest/gtest.h"
#include "iostream"
#include "Model.h"
#include "Result.h"
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>

/*
 * Test for performing Breusch-Pagan model
 * */

static std::vector<double> DOSAGES = {
    0, 2, 0, 2, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 2, 2,
    0, 1, 0, 1, 2, 0, 1, 1, 1, 1, 1, 0, 2, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 2, 0, 2, 1, 2, 1, 0, 2, 0, 0, 1, 1, 1, 0,
    2, 0, 0, 1, 2, 0, 1, 1, 1, 1, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 2, 0, 0, 2, 0, 2, 1, 1, 0, 1, 0, 1,
    0, 2, 2, 0, 1, 1, 0, 0, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 0, 2, 2, 0, 1, 1, 0, 0, 0, 1, 0, 1, 2, 2, 0, 2, 1, 0, 0,
    0, 1, 1, 0, 0, 1, 2, 0, 0, 0, 2, 0, 1, 0, 2, 1, 0, 1, 2, 1, 0, 1, 1, 0, 0, 2, 0, 0, 1, 1, 1, 1, 2, 0, 1, 2, 1, 1,
    0, 2, 1, 1, 2, 1, 1, 0, 1, 1
};
static std::vector<double> PHENO = {
    0.522282168524833, -5.97159895922116, -0.440526196202908, -5.73629947937107, 0.565161742547473, -1.19646012034654,
    0.595643003888263, 0.343958352102906, -0.329072974464988, 4.05568050968734, -4.11814091986007, -3.5537295154718,
    1.32029371983474, 5.7493407446655, 2.21709180256317, -1.46274426579031, -0.439619003751496, 0.657412424396651,
    8.72368664203702, -0.595398935789, 5.03774698616385, 4.62997815935095, 3.41031850809713, 0.474943044148217,
    -2.18594644011259, -2.0648595591344, 0.475435421232725, 0.390281600028403, -6.45700963893647, 0.986555948076592,
    1.42398359867007, 3.26994917129858, 2.11857253140262, 0.860124249616191, 8.20762563718038, -2.31697274908362,
    -6.67888590076259, 5.90420665016735, 0.499881065667666, -0.452503547529653, 2.07670993327744, -3.68765540931014,
    -1.8491082440702, 0.532788293060056, 0.232497686539817, 1.6705767524884, -5.28103291566824, 3.23426068062678,
    0.238355888856576, 0.77890003053602, -8.31101598685559, -0.0867079053484027, -0.630268076943675,
    -3.08362852090857, 1.05773736712652, -0.360320800556423, 0.350593765015391, -0.663474385430901, 0.94809485288954,
    -2.47809290925924, -0.375823436339545, -0.671311036017151, -0.515398662055499, -0.189622581277104,
    -5.03241940707449, 0.138922205901387, 4.46133220551349, -1.55903846958601, 0.672614533611014, 0.154682805363378,
    -0.146026391952837, 1.70143183167595, 0.856649961275889, 3.19965154977335, -3.87194465181213, -0.77890222515892,
    15.7939690036799, -1.96454419520035, 0.769170668260007, -2.5982675338341, -5.0741526508138, -0.102580508795901,
    -0.496008464068098, 2.94839306662081, 3.09874033867321, 0.34896910759503, 0.904243556201491, 2.24203613746722,
    -1.19512289061954, -1.55738814284606, 4.44251709196676, 0.49817216339198, 2.63052250106299, -0.367202719244542,
    0.262331084814879, 0.262704124824336, 0.640738877993306, 0.307089830771924, -1.01502221189327, -1.05392730103636,
    -3.15533145405574, 0.00214395093177762, 2.2297414814471, -1.00177908621909, -0.61722192863403, 1.99030875320093,
    -0.084818645674312, 3.78561380875635, -0.646734676930393, -0.954808336004956, -0.971229202960724,
    -2.31073189265367, 1.20468703481302, 6.20047536362092, -0.882457683653546, 2.33369568809647, 4.63355788045797,
    -0.930628616693406, -6.89245862400824, 1.56459985012352, 0.514242102223054, 0.293168380759209, -0.704768631523196,
    0.207142723762983, 1.2280774574073, -0.430012569304131, 0.712936016878738, -0.445622389981341, 5.85932073306511,
    0.609058662925924, -0.317294209895932, -0.0476483564726301, 0.170678010698749, -1.22507132610198, 5.9981240427872,
    -1.59268405931497, 2.23912608412108, -0.641714188063126, 2.85933535125804, 1.93544276888011, 0.394249899631979,
    0.632521494371694, -0.354395798676209, 1.01021099822561, 1.19664554149654, 5.43393113609091, 15.2153783088347,
    -2.19720851198184, 4.08306883492993, 1.38546168757324, 0.713078103934547, 1.02826853911963, -0.570632725418573,
    -4.21193842733121, 0.537256877476888, 1.20354481603093, 1.49618269695644, -1.17473695943517, -4.85762413074619,
    0.352538646178498, -1.28618864655408, -0.272522855069837, -0.0874495918047668, -0.88181286823508,
    4.36005856300851, -0.855250672623498, -0.203229680356456, 0.736925075090769, -2.07823365938166, -1.73640616601282,
    -1.72808470104198, -1.73382710673325, -2.23662183763668, -5.15500700324523, 1.85275781893222, -0.497579143100874,
    -0.601027404256258, -1.71127715230004, -0.00493237975734888, -0.0971479127091374, -4.02739153888159,
    0.620921307289942, -2.31459048673085, 5.24840330286662, 2.82164593721117, 0.210930786953303, 1.32541204793777,
    -3.11935648433972, 3.19528875708452, 2.65183102290013, -1.16278371681341, -2.62199741751196, -0.498018360244571,
    -0.600628207986387, -5.1989262441418, 1.02757089200525, -6.54227918659197, -0.846710740540237, -1.99565028281955,
    1.28624046292148
};

TEST(ModelTest, fit) {
  std::vector<double> pheno = PHENO;
  std::vector<double> dosages = DOSAGES;
  assert(dosages.size() == pheno.size());
  int n = pheno.size();
  int p = 1; // n of covariates
  Eigen::MatrixXd X1 = Eigen::MatrixXd(n, p + 1); // add 1 for intercept
  Eigen::MatrixXd X2 = Eigen::MatrixXd(n, p * 2 + 1); // added xsq and intercept
  Eigen::VectorXd y = Eigen::VectorXd(n);

  // initialise empty matrix
  for (unsigned i = 0; i < n; ++i) {
    X1(i, 0) = 1; // intercept
    X1(i, 1) = 0; // X
    X2(i, 0) = 1; // intercept
    X2(i, 1) = 0; // X
    X2(i, 2) = 0; // Xsq
    y(i, 0) = pheno[i];
  }

  // fit B-P model
  std::string chr = "01";
  std::string rsid = "rs1";
  std::string allele = "A";
  std::set<unsigned> non_nulls_idx;
  for (unsigned i = 0; i < dosages.size(); i++) {
    non_nulls_idx.insert(i);
  }
  jlst::Result result = jlst::Model::fit(chr, 1, rsid, allele, allele, dosages, non_nulls_idx, X1, X2, y);

  // check estimate and SE are similar to R
  ASSERT_NEAR(result.beta, 0.262569, 0.01);
  ASSERT_NEAR(result.se, 0.301703, 0.01);
  ASSERT_NEAR(result.pval, 0.385, 0.01);
  ASSERT_NEAR(result.t, 0.870, 0.01);
  ASSERT_NEAR(result.phi_x, 3.743, 0.01);
  ASSERT_NEAR(result.se_x, 6.797, 0.01);
  ASSERT_NEAR(result.phi_xsq, 4.956, 0.01);
  ASSERT_NEAR(result.se_xsq, 3.504, 0.01);
  ASSERT_NEAR(result.phi_pval, 2.225e-07, 0.01);
  ASSERT_NEAR(result.phi_f, 16.573, 0.01);
}

TEST(ModelTest, fit_missing_vals) {
  std::vector<double> pheno = PHENO;
  std::vector<double> dosages = DOSAGES;
  assert(dosages.size() == pheno.size());
  dosages[0] = -1;
  int n = pheno.size();
  int p = 1;
  Eigen::MatrixXd X1 = Eigen::MatrixXd(n, p + 1);
  Eigen::MatrixXd X2 = Eigen::MatrixXd(n, p * 2 + 1);
  Eigen::VectorXd y = Eigen::VectorXd(n);

  // initialise empty matrix
  for (unsigned i = 0; i < n; ++i) {
    X1(i, 0) = 1; // intercept
    X1(i, 1) = 0; // X
    X2(i, 0) = 1; // intercept
    X2(i, 1) = 0; // X
    X2(i, 2) = 0; // Xsq
    y(i, 0) = pheno[i];
  }

  // fit B-P model
  std::string chr = "01";
  std::string rsid = "rs1";
  std::string allele = "A";
  std::set<unsigned> non_nulls_idx;
  for (unsigned i = 0; i < dosages.size(); i++) {
    non_nulls_idx.insert(i);
  }
  non_nulls_idx.erase(1);
  jlst::Result result = jlst::Model::fit(chr, 1, rsid, allele, allele, dosages, non_nulls_idx, X1, X2, y);
  ASSERT_EQ(result.n, dosages.size() - 2);
}

TEST(ModelTest, ftest) {
  unsigned n = 321305;
  unsigned df_f = n - 3;
  unsigned df_r = n - 1;
  unsigned df_n = df_r - df_f;
  double f = 500;
  boost::math::fisher_f dist(df_n, df_f);
  double pval = boost::math::cdf(boost::math::complement(dist, f));
  //printf("%e\n", pval);
  ASSERT_NEAR(pval, 1.548735e-217, 1.548735e-217 * .001);
}

TEST(ModelTest, ttest) {
  int df = 198;
  double t = 0.87;
  double p = 0.385;
  boost::math::students_t t_dist(df);
  double pval = 2.0 * boost::math::cdf(boost::math::complement(t_dist, fabs(t)));
  ASSERT_NEAR(pval, p, p * .001);
}