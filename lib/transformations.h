#ifndef TRANSFORMATION_LIB_TRANSFORMATIONS_H_
#define TRANSFORMATION_LIB_TRANSFORMATIONS_H_

#include "radian_degree.h"

#include <array>
#include <string>

enum class ELLIPSOID { PZ90, WGS84, SK42 };

struct Params {
    double a;
    double e2;
    double da;
    double de2;
    double dx;
    double dy;
    double dz;
};

class Geo {
 protected:
    static double dB(Radian B, Radian L, double H, Params p);
    static double dL(Radian B, Radian L, double H, Params p);
    static double dH(Radian B, Radian L, double H, Params p);

    static constexpr double ro = 206264.8062;
};

class WGS84;
class SK42;
class PZ90;
class GaussKruger;
class UTM;

class WGS84 : public Geo {
 public:
    WGS84() = default;
    WGS84(Degree latitude, Degree longitude, double altitude)
        : latitude(latitude), longitude(longitude), altitude(altitude) {}
    explicit WGS84(SK42 sk_42);
    explicit WGS84(PZ90 pz_90);
    explicit WGS84(UTM utm);
    Degree latitude{};
    Degree longitude{};
    double altitude{};

    static constexpr double _a = 6378245;
    static constexpr double _al = 1 / 298.3;
    static constexpr double _e2 = 2 * _al - _al * _al;
};

class SK42 : public Geo {
 public:
    explicit SK42(WGS84 wgs_84);
    explicit SK42(GaussKruger gk);

    Degree latitude{};
    Degree longitude{};
    double altitude{};
    Params p {
        (_a + WGS84::_a) / 2,
        (_e2 + WGS84::_e2) / 2,
        WGS84::_a - _a,
        WGS84::_e2 - _e2,
        23.92,
        -141.27,
        -80.9
    };

 private:
    static constexpr double _a = 6378137;
    static constexpr double _al = 1 / 298.257223563;
    static constexpr double _e2 = 2 * _al - _al * _al;
};

class PZ90 : public Geo {
 public:
    PZ90(Degree latitude, Degree longitude, double altitude);
    explicit PZ90(WGS84 wgs_84);

    Degree latitude{};
    Degree longitude{};
    double altitude{};
    Params p {
        (_a + WGS84::_a) / 2,
        (_e2 + WGS84::_e2) / 2,
        WGS84::_a - _a,
        WGS84::_e2 - _e2,
        -1.1,
        -0.3,
        -0.9
    };

 private:
    static constexpr double _a = 6378136.5;
    static constexpr double _al = 1 / 298.25784;
    static constexpr double _e2 = 2 * _al - _al * _al;
};

class UTM {
 public:
    explicit UTM(WGS84 wgs_84);
    UTM(Degree E, Degree N, double altitude, std::string  zone);

    double E{};
    double N{};
    double altitude{};
    std::string zone{};

    static constexpr double k0 = 0.9996;
    static constexpr double E0 = 500000.0;
    static constexpr double N0 = 10000000.0;

 private:
    static char letter_designator(Degree latitude);
};

class GaussKruger {
 public:
    GaussKruger() = default;
    explicit GaussKruger(SK42 sk_42);

    double x;
    double y;
    double height;
};

#endif  // TRANSFORMATION_LIB_TRANSFORMATIONS_H_
