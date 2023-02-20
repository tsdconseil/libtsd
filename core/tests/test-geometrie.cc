#include "tsd/tsd.hpp"
#include "tsd/geometrie.hpp"
#include "tsd/vue.hpp"
#include "tsd/tests.hpp"


using namespace tsd;
using namespace tsd::geo;
using namespace tsd::vue;
using namespace tsd::temps;

using Eigen::Matrix3f;
using Eigen::Vector3f;
using Eigen::Vector2f;

template<typename T1, typename T2>
  bouléen sim(T1 a, T2 b)
{
  retourne abs(a - b) < 1e-7;
}





void test_quaternions()
{
  msg("Test quaternions...");

  {
    Quaternion q = Quaternion::identite();
    soit m = q.rot_mat();
    tsd_assert((m - Eigen::Matrix3f::Identity()).norm() < 1e-6);
    tsd_assert((q.inv().rot_mat() - Eigen::Matrix3f::Identity()).norm() < 1e-6);
  }

  {
    Eigen::Array4f a = Eigen::Array4f::Random();
    Quaternion q1(a);
    soit R = q1.rot_mat();

    pour(auto k = 0; k < 3; k++)
    {
      soit n = R.col(k).norm();
      msg("R: norme colonne {} = {}", k, n);
      tsd_assert_msg(n - 1 < 1e-5, "Matrice de rotation invalide.");
    }

    Quaternion q2(R);

    msg("q1 = {}, q2 = {}", q1.q, q2.q);
    msg("q1.norm() = {}, q2.norm() = {}", q1.q.abs2().sum(), q2.q.abs2().sum());
    soit err = sqrt((q1.q - q2.q).abs2().mean());
    err = min(err, sqrt((q1.q + q2.q).abs2().mean()));
    tsd_assert_msg(err < 1e-5, "Erreur trop importante ({})", err);


    Eigen::Vector3f x  = Eigen::Vector3f::Random();
    Eigen::Vector3f y1 = R * x;
    Eigen::Vector3f y2 = q1.rotate(x);
    Eigen::Vector3f y3 = q2.rotate(x);

    msg("x  = {}", x);
    msg("y1 = {}", y1);
    msg("y2 = {}", y2);
    msg("y3 = {}", y3);

    tsd_assert_msg(y1.isApprox(y2, 1e-5) && y1.isApprox(y3, 1e-5), "Pb rotation quat.");
  }

  msg("ok");
}



entier test_geometrie()
{
  test_quaternions();
  retourne 0;
}
