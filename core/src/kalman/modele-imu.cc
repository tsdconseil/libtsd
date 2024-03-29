#include "tsd/tsd.hpp"
#include "tsd/apps/kalman.hpp"


using namespace tsd::geo;


namespace tsd::kalman {

// IMU (Inertial Measuring Unit) generic model
//
// Calling Sequence
// ssm = model_imu();
// ssm = model_imu(cfg);
//
// Parameters
// cfg: IMU model configuration (see <link linkend="model_imu_cfg">model_imu_cfg()</link>)
// ssm: State / space model
//
// Description
// Generic IMU model with Euler angles parametrization, and with the following
// (optionnal) sensors:
// <itemizedlist>
//   <listitem><para>Accelerometer (3 axis)</para></listitem>
//   <listitem><para>Gyroscope (3 axis)</para></listitem>
//   <listitem><para>Magnetometer / compas (3 axis)</para></listitem>
//   <listitem><para>Barometer</para></listitem>
//   <listitem><para>GPS</para></listitem>
// </itemizedlist>
//
// This model can be used either pour simulation (see <link linkend="ssm_simu">ssm_simu</link>)
// or pour orientation / position estimation from experimental data (pour instance using an EKF filter, see <link linkend="ekf_init">ekf_init</link>).
//
//  The model can use as a state vector either only the object orientation, or both the
//  object orientation and localization.
// So it is two parts, where each part is optionnal (according to the configuration):
// <itemizedlist>
//   <listitem><para><emphasis role="bold">Position: </emphasis>  <latex>$(x,y,z,\dot{x},\dot{y},\dot{z},\ddot{x},\ddot{y},\ddot{z})$</latex></para></listitem>
//   <listitem><para><emphasis role="bold">Orientation: </emphasis>   <latex>$(\phi,\theta,\psi,\dot{\phi},\dot{\theta},\dot{\psi})$</latex></para></listitem>
// </itemizedlist>
//
// The Observation vector is in this order: <latex>$(a_x,a_y,a_z,\omega_x,\omega_y,\omega_y,m_x,m_y,m_z,p_{baro},x_{gps},y_{gps},z_{gps})$</latex>,
//  where each element is optionnal, according to the configuration.
//


  // Configuration pour IMU model
  //
  // Calling Sequence
  // cfg = model_imu_cfg(); // Default configuration
  //
  // Parameters
  // cfg: IMU configuration structure
  //
  // Description
  // The configuration structure fields can be changed afterwards. Available fields are:
  // <itemizedlist>
  //   <listitem><para><emphasis role="bold">Specification of states to compute:</emphasis></para></listitem>
  // <itemizedlist>
  //   <listitem><para><varname>cfg.estim_orient</varname>: If oui, the model includes orientation (Euler angles) estimation. </para></listitem>
  //   <listitem><para><varname>cfg.estim_pos</varname>: If oui, the model includes position estimation. </para></listitem>
  // </itemizedlist>
  //   <listitem><para><emphasis role="bold">Specification of available sensors:</emphasis></para></listitem>
  // <itemizedlist>
  //   <listitem><para><varname>cfg.acc_available</varname>: Boolean specifying si accelerometer measurements are available or not</para></listitem>
  //   <listitem><para><varname>cfg.gyro_available</varname>: Boolean specifying si gyroscope measurements are available or not</para></listitem>
  //   <listitem><para><varname>cfg.mag_available</varname>: Boolean specifying si magnetometer (compas) measurements are available or not</para></listitem>
  //   <listitem><para><varname>cfg.baro_available</varname>: Boolean specifying si barometer measurements are available or not</para></listitem>
  //   <listitem><para><varname>cfg.gps_available</varname>: Boolean specifying si GPS measurements are available or not</para></listitem>
  // </itemizedlist>
  //   <listitem><para><emphasis role="bold">Specification of sensors noise:</emphasis></para></listitem>
  // <itemizedlist>
  //   <listitem><para><varname>cfg.anoise</varname>: Noise standard deviation on the accelerometers outputs (scalar)</para></listitem>
  //   <listitem><para><varname>cfg.gnoise</varname>: Noise standard deviation on the gyroscopes outputs (scalar)</para></listitem>
  //   <listitem><para><varname>cfg.mnoise</varname>: Noise standard deviation on the compas outputs (scalar)</para></listitem>
  //   <listitem><para><varname>cfg.bnoise</varname>: Noise standard deviation on the barometer output (scalar)</para></listitem>
  //   <listitem><para><varname>cfg.gpsnoise</varname>: Noise standard deviation on the GPS outputs (scalar)</para></listitem>
  // </itemizedlist>
  // </itemizedlist>
  //
  // <refsection><title>Example 1: orientation estimation from accelerometer and gyroscopes</title></refsection>
  // <programlisting>
  //   cfg = model_imu_cfg();
  //   cfg.estim_orient   = %t;
  //   cfg.acc_available  = %t;
  //   cfg.gyro_available = %t;
  //   sys = model_imu(cfg);
  //   // sys can now be used with ekf_init (pour estimation)
  //   // or ssm_simu (pour simulation)
  // </programlisting>
  //
  // <refsection><title>Example 2: orientation and position estimation from accelerometer, gyroscopes and GPS</title></refsection>
  // <programlisting>
  //   cfg = model_imu_cfg();
  //   cfg.estim_orient   = %t;
  //   cfg.estim_pos      = %t;
  //   cfg.acc_available  = %t;
  //   cfg.gyro_available = %t;
  //   cfg.gps_available  = %t;
  //   sys = model_imu(cfg);
  //   // sys can now be used with ekf_init (pour estimation)
  //   // or ssm_simu (pour simulation)
  // </programlisting>
  //








struct IMU: SSM
{
  IMUConfig config;

  float G  = 9.81;
  entier ipos, iorient, iacc, igyro, imag, ibaro, igps;
  float δt;

  IMU(const IMUConfig &cfg)
  {
    ns = 0;
    config = cfg;
    δt = 1.0 / cfg.fs;
    si(cfg.estim_pos)
    {
      ipos = ns;
      ns += 9; // 3 pos + dérivées + acc
    }
    si(cfg.estim_orient)
    {
      iorient = ns;
      ns += 6; // 3 angles + dérivées
    }


    no = 0;
    si(cfg.acc_available)
    {
      iacc = no;
      no += 3;
    }
    si(cfg.gyro_available)
    {
      igyro = no;
      no += 3;
    }
    si(cfg.mag_available)
    {
      imag = no;
      no += 3;
    }
    si(cfg.baro_available)
    {
      ibaro = no;
      no += 1;
    }
    si(cfg.gps_available)
    {
      igps = no;
      no += 3;
    }

    // Bruit de process :
    //  = 0 sur x,y,z,phi,theta,psi,
    //  = 1 sur dphi, ...
    //  = 0 sur dx,dy,...
    //  = 10m/s^-2 sur les accélérations linéaires
    //Q = diag([1e-10 * ones(6,1) ; ones(3,1) ; 1e-10 * ones(3,1) ; 10^2*ones(3,1)]);

    Vecf v(ns);

    si(cfg.estim_pos)
    {
      v.head(6).setConstant(1e-5f);
      v.segment(6,3).setConstant(10.0f);
    }

    v.tail(3).setConstant(1.0f);
    v.segment(v.rows() - 6, 3).setConstant(1e-5f);

    //Q = diag([v ; 1e-5 * ones(3,1) ; 1*ones(3,1)])^2;

    //v = [];

    v.resize(no);


    // Bruit de mesure
    si(cfg.acc_available)
      v.head(3).setConstant(cfg.anoise);
    si(cfg.gyro_available)
      v.segment(igyro, 3).setConstant(cfg.gnoise);
    si(cfg.mag_available)
      v.segment(imag, 3).setConstant(cfg.mnoise);
    si(cfg.baro_available)
      v(ibaro) = cfg.bnoise;
    si(cfg.gps_available)
      v.segment(igps, 3).setConstant(cfg.gpsnoise);

    R = Tabf::diagonal(abs2(v));

    //ssm = ssm_init(ns,no,f,g,[],[],Q,R);
    //ssm.cfg = cfg;

    noms_etats.clear();
    si(cfg.estim_pos)
      noms_etats = {"x","y","z","dx/dt","dy/dt","dz/dt","d^2x/dt^2","d^2y/dt^2","d^2z/dt^2"};

    noms_etats.push_back("ϕ");
    noms_etats.push_back("θ");
    noms_etats.push_back("ψ");
    noms_etats.push_back("dϕ/dt");
    noms_etats.push_back("dθ/dt");
    noms_etats.push_back("dψ/dt");

    noms_obs.clear();

    si(cfg.acc_available)
    {
      noms_obs.push_back("ax");
      noms_obs.push_back("ay");
      noms_obs.push_back("az");
    }
    si(cfg.gyro_available)
    {
      noms_obs.push_back("ωx");
      noms_obs.push_back("ωy");
      noms_obs.push_back("ωz");
    }
    si(cfg.mag_available)
    {
      noms_obs.push_back("mx");
      noms_obs.push_back("my");
      noms_obs.push_back("mz");
    }
    si(cfg.baro_available)
      noms_obs.push_back("zbaro");
    si(cfg.gps_available)
    {
      noms_obs.push_back("xgps");
      noms_obs.push_back("ygps");
      noms_obs.push_back("zgps");
    }


    //ssm2 = ssm_init(ssm.ns,ssm.no,f,ssm.g,[],[],Q,ssm.R);
    //ssm2.cfg = ssm.cfg;

    // Same observation noise
    //R = zeros(ssm.ns,ssm.ns);

    // Modèle IMU avec quaternion / sans quaternion
    // Ici modèle IMU avec angles d'Euler
  }

  Vecf f(const Vecf &x) const
  {
    Vecf x2(ns);

    si(config.estim_pos)
    {
      x2.head(3) = x.head(3) + δt * x.segment(3, 3);
      x2.segment(3,3) = x.segment(3,3) + δt * x.segment(6,3);
      x2.segment(6,3).setZero();//x(7:9); // Accelerations nulles (bruit sans mémoire)
    }

    si(config.estim_orient)
    {
      entier io = iorient;
      x2.segment(io,3) = x.segment(io,3) + δt * x.segment(io+3,3);
      x2.segment(io+3,3) = x.segment(io+3,3); // Vitesses angulaires constantes
    }
    retourne x2;
  }

  Vecf g(const Vecf &x) const
  {
    soit io = iorient, id = 0;
    soit cp = cos(x(io)),
         sp = sin(x(io)),
         ct = cos(x(io+1)),
         st = sin(x(io+1));

    // Matrice de rotation repère absolu -> repère inertiel
    soit Rmat = (Cardan{x(io),x(io+1),x(io+2)}).mat_rotation();

    Vecf y(no);

    si(config.acc_available)
    {
      soit v = Vecf::zeros(3);
      v(2) = G;

      si(config.estim_pos)
        v += x.segment(6, 3);

      y.head(3) = Rmat.matprod(v);
      id += 3;
    }
    si(config.gyro_available)
    {
      soit dϕ = x(io+3), dθ = x(io+4), dψ = x(io+5);
      y(id)   = dϕ - st * dψ;
      y(id+1) = cp * dθ + sp * ct * dψ;
      y(id+2) = cp * ct * dψ - sp * dθ;
      id += 3;
    }
    si(config.mag_available)
    {
      soit v = Vecf::valeurs({0, 1, 0});
      y.segment(id,3) = Rmat.matprod(v); // En supposant force normalisée à 1 gauss.
      id += 3;
    }
    si(config.baro_available)
    {
      y(id) = x(2);
      id++;
    }
    si(config.gps_available)
    {
      y.segment(id,3) = x.head(3);
      id++;
    }
    retourne y;
  }







#if 0
function ssm2 = model_imu_for_simu(ssm)

    //ssm2 = model_imu(ssm.cfg);

    // Process noise is null
    Q = zeros(ssm.ns,ssm.ns);

function x2 = f(x,ssm)
    δt = ssm.cfg.δt;

    si(cfg.estim_pos)
        x2(1:3) = x(1:3) + δt * x(4:6);
        x2(4:6) = x(4:6) + δt * x(7:9);
        x2(7:9) = 0;//x(7:9); // Accelerations nulles (bruit sans mémoire)
    end;

    si(cfg.estim_orient)
        io = ssm.cfg.iorient;
        x2(io:io+2) = x(io:io+2) + δt * x(io+3:io+5);
        x2(io+3:io+5) = x(io+3:io+5); // Vitesses angulaires constantes

        x2(io+3) = %pi/2 + ssm.k * %pi/30;
        x2(io+4) = %pi/4 + ssm.k * %pi/10;
        x2(io+5) = -%pi/4 - ssm.k * %pi/30;
    end;



    //x2(7) = 10;

endfunction;
#endif



};


sptr<SSM> modele_imu(const IMUConfig &cfg)
{
  retourne std::make_shared<IMU>(cfg);
}

}



