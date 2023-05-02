#include "tsd/tsd.hpp"
#include "tsd/geometrie.hpp"
#include <cmath>
#include <iomanip>
#include "Eigen/Geometry"


using namespace std;
using namespace tsd::temps;


ostream_formater(Eigen::Array4f)

namespace tsd::geo {




ostream& operator<<(ostream& ss, const Cardan &t)
{
  ss << sformat("Cardan[φ={},θ={},ψ={}]", t.φ, t.θ, t.ψ);
  retourne ss;
}

ostream& operator<<(ostream& ss, const Quaternion &t)
{
  ss << sformat("Quaternion[{}]", t.q);
  retourne ss;
}

static Eigen::Matrix3f makeC(const Eigen::Vector3f &v)
{
  Eigen::Matrix3f R(3,3);
  R(0,0) = 0;
  R(0,1) = -v(2);
  R(0,2) = v(1);

  R(1,0) = v(2);
  R(1,1) = 0;
  R(1,2) = -v(0);

  R(2,0) = -v(1);
  R(2,1) = v(0);
  R(2,2) = 0;

  retourne R;
}



// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
Eigen::Matrix3f Quaternion::rot_mat() const
{
  Eigen::Matrix3f R;

  soit w = q(0), x = q(1), y = q(2), z = q(3);

  R << 1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w),
       2*(x*y+z*w), 1-2*(z*z+x*x), 2*(y*z-x*w),
       2*(x*z-y*w), 2*(y*z+x*w), 1-2*(y*y+x*x);

  retourne R;
}

static inline Eigen::Vector3f cross(const Eigen::Vector3f &a, const Eigen::Vector3f &b)
{
  retourne Eigen::Vector3f(
      a(1) * b(2) - a(2) * b(1),
      a(2) * b(0) - a(0) * b(2),
      a(0) * b(1) - a(1) * b(0));
}

// W = V + ow_fastcross(2*Q(1:3,:),ow_fastcross(Q(1:3,:),V) - Q(4,:).*V);

Eigen::Vector3f Quaternion::rotate(const Eigen::Vector3f &x) const
{
  Eigen::Vector3f e = -q.tail(3).matrix();
  float s = q(0);
  retourne x + cross(2 * e, cross(e, x) - s * x);
  //retourne rot_mat() * x;
}

Eigen::Matrix4f Quaternion::mat() const
{
  Eigen::Matrix4f R;

  R(0,0) = q(0);
  R.block(0, 1, 1, 3) = -q.tail(3).transpose();
  R.block(1, 0, 3, 1) = q.tail(3);
  R.block(1, 1, 3, 3) = makeC(q.tail(3));

  retourne R;
}


Quaternion::Quaternion(const Eigen::Array4f &q_)
{
  q = q_ / sqrt(q_.abs2().sum());
}

Quaternion::Quaternion(float q0, float q1, float q2, float q3)
{
  q << q0, q1, q2, q3;
  q /= sqrt(q.abs2().sum());
}


Quaternion::Quaternion(const Eigen::Matrix3f &R)
{
  // D'après https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
  float t;

  si(R(2,2) < 0)
  {
    si(R(0,0) > R(1,1))
    {
      t = 1 + R(0,0) - R(1,1) - R(2,2);
      q << R(2,1)-R(1,2), t, R(1,0)+R(0,1), R(2,0)+R(0,2);
    }
    sinon
    {
      t = 1 - R(0,0) + R(1,1) - R(2,2);
      q << R(0,2)-R(2,0), R(0,1) + R(1,0), t, R(1,2)+R(2,1);
    }
  }
  sinon
  {
    si(R(0,0) < -R(1,1))
    {
      t = 1 - R(0,0) - R(1,1) + R(2,2);
      q << R(1,0)-R(0,1), R(2,0)+R(0,2), R(1,2)+R(2,1), t;
    }
    sinon
    {
      t = 1 + R(0,0) + R(1,1) + R(2,2);
      q << t, R(2,1)-R(1,2), R(0,2)-R(2,0), R(1,0)-R(0,1);
    }
  }

  q *= 0.5 / sqrt(t);
}


Quaternion Quaternion::identite()
{
  retourne Quaternion{1,0,0,0};
}


Quaternion Quaternion::inv() const
{
  Quaternion res;
  res.q(0) = q(0);
  res.q.tail(3) = - q.tail(3);
  retourne res;
}


Cardan::Cardan(float φ, float θ, float ψ)
{
  this->φ = φ;
  this->θ = θ;
  this->ψ = ψ;
}

Cardan::Cardan(const Eigen::Matrix3f &R)
{
  φ = atan2(R(1,2),R(2,2));
  θ = -asin(R(0,2));
  ψ = atan2(R(0,1),R(0,0));
}

Cardan::Cardan(const Quaternion &q)
{
  soit v = q.q;
  // Conversion from quaternion
  // eqn (290) from "Representing Attitude: Euler Angles, Unit Quaternions,
  // and Rotation Vectors", James Diebel, Oct 2006
  φ = atan2(2*v(2)*v(3)+2*v(0)*v(1), v(3)*v(3)-v(2)*v(2)-v(1)*v(1)+v(0)*v(0));
  θ = -asin(2*v(1)*v(3)-2*v(0)*v(2));
  ψ = atan2(2*v(1)*v(2)+2*v(0)*v(3), v(1)*v(1)+v(0)*v(0)-v(3)*v(3)-v(2)*v(2));
}

/*static Eigen::Matrix3f rotmat_3d_R1(float α)
{
  float ca = cos(α), sa = sin(α);
  Eigen::Matrix3f R;
  R << 1, 0, 0,
       0, ca, sa,
       0, -sa, ca;
  retourne R;
}

static Eigen::Matrix3f rotmat_3d_R2(float α)
{
  float ca = cos(α), sa = sin(α);
  Eigen::Matrix3f R;
  R << ca, 0, -sa,
       0, 1, 0,
       sa, 0, ca;
  retourne R;
}*/




/*static Eigen::Matrix3f rotmat_3d_R3(float α)
{
  float ca = cos(α), sa = sin(α);
  Eigen::Matrix3f R;
  R << ca, sa, 0,
       -sa, ca, 0,
       0, 0, 1;
  retourne R;
}*/

/*Eigen::Matrix3f rotmat_3d(float α, entier axe)
{
  si(axe == 0)
    retourne rotmat_3d_R1(α);
  sinon si(axe == 1)
    retourne rotmat_3d_R2(α);
  retourne rotmat_3d_R3(α);
}*/


Tabf Cardan::mat_rotation() const
{
  Eigen::Matrix3f r = rotmat_3d_R1(φ) * rotmat_3d_R2(θ) * rotmat_3d_R3(ψ);
  retourne (Tabf::map(r.data(), 3, 3)).clone();
}


// Examples
// c = cardan(π/4,0,0) // a rotation of angle phi = π/4 (roll only)
// c.phi                 // extract roll angle
// c.theta               // extract pitch angle
// c.psi                 // extract yaw angle



// Examples
//  // R = rotation matrix, pour phi = %pi/4 and psi = %pi/2
//  R = rotmat(cardan(%pi/4,0,%pi/2));




}

