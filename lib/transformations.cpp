#include "transformations.h"

#include <cmath>
#include <utility>

using std::sin;
using std::cos;
using std::sinh;
using std::cosh;
using std::tanh;
using std::atanh;
using std::pow;
using std::round;
using std::sqrt;

double Geo::dB(Radian B, Radian L, double H, Params p) {
    double M = p.a * (1 - p.e2) / pow((1 - p.e2 * pow(sin(B), 2)), 1.5);
    double N = p.a * pow((1 - p.e2 * pow(sin(B), 2)), -0.5);
    return ro / (M + H) * (N / p.a * p.e2 * sin(B) * cos(B) * p.da + ((N * N) / (p.a * p.a) + 1) * N * sin(B) * cos(B) * p.de2 / 2 - (p.dx * cos(L) + p.dy * sin(L)) * sin(B) + p.dz * cos(B));
}
double Geo::dL(Radian B, Radian L, double H, Params p) {
    double N = p.a * pow((1 - p.e2 * pow(sin(B), 2)), -0.5);
    return ro / ((N + H) * cos(B)) * (-p.dx * sin(L) + p.dy * cos(L));
}
double Geo::dH(Radian B, Radian L, double H, Params p) {
    double N = p.a * pow((1 - p.e2 * pow(sin(B), 2)), -0.5);
    double dH = -p.a / N * p.da + N * pow(sin(B), 2) * p.de2 / 2 + (p.dx * cos(L) + p.dy * sin(L)) * cos(B) + p.dz * sin(B);
    return dH;
}

WGS84::WGS84(SK42 sk_42) {
    altitude = sk_42.altitude;
    latitude = Degree{sk_42.latitude + dB(sk_42.latitude, sk_42.longitude, sk_42.altitude, sk_42.p) / 3600};
    longitude = Degree{sk_42.longitude + dL(sk_42.latitude, sk_42.longitude, sk_42.altitude, sk_42.p) / 3600};
}
WGS84::WGS84(PZ90 pz_90) {
    altitude = pz_90.altitude;
    latitude = Degree{pz_90.latitude + dB(pz_90.latitude, pz_90.longitude, pz_90.altitude, pz_90.p) / 3600};
    longitude = Degree{pz_90.longitude + dL(pz_90.latitude, pz_90.longitude, pz_90.altitude, pz_90.p) / 3600};
}
WGS84::WGS84(UTM utm) {
    altitude = utm.altitude;
    double e1 = (1 - sqrt(1 - WGS84::_e2)) / (1 + sqrt(1 - WGS84::_e2));

    double x = utm.E - UTM::E0;  // remove 500,000 meter offset for longitude
    double y = utm.N;

    char zoneLetter = utm.zone[utm.zone.length() - 1];
    int zoneNumber = std::stoi(utm.zone.substr(0, utm.zone.length() - 1));
    if ((zoneLetter - 'N') < 0) {
        // remove 10,000,000 meter offset used for southern hemisphere
        y -= UTM::N0;
    }

    // +3 puts origin in middle of zone
    double lambda0 = (zoneNumber - 1) * 6 - 180 + 3;
    static const double EPrimeSquared = (WGS84::_e2) / (1 - WGS84::_e2);

    double M = y / UTM::k0;
    double mu = M /
        (WGS84::_a * (1 - WGS84::_e2 / 4 - 3 * WGS84::_e2 * WGS84::_e2 / 64
            - 5 * WGS84::_e2 * WGS84::_e2 * WGS84::_e2 / 256));

    Radian phi1Rad = Radian{mu + sin(2 * mu) * (3 * e1 / 2 - (27. / 32) * pow(e1, 3) / 32) +
        sin(4 * mu) * ((21. / 16) * pow(e1, 2) - (55. / 32) * pow(e1, 4)) +
        sin(6 * mu) * ((151. / 96) * pow(e1, 3))};

    double N1 = WGS84::_a / sqrt(1 - WGS84::_e2 * pow(sin(phi1Rad), 2));
    double T1 = tan(phi1Rad) * tan(phi1Rad);
    double C1 = EPrimeSquared * pow(cos(phi1Rad), 2);
    double R1 = WGS84::_a * (1 - WGS84::_e2) / pow(1 - WGS84::_e2 * pow(sin(phi1Rad), 2), 1.5);
    double D = x / (N1 * UTM::k0);

    Radian latRad = Radian{phi1Rad - ((N1 * tan(phi1Rad) / R1) *
        (D * D / 2 - (5 + 3 * T1 + 10 * C1 - 4 * pow(C1, 2) - 9 * EPrimeSquared) * pow(D, 4) / 24 +
            (61 + 90 * T1 + 298 * C1 + 45 * pow(T1, 2) - 252 * EPrimeSquared - 3 * pow(C1, 2))
                * pow(D, 6) / 720))};

    latitude = latRad;

    Radian longRad = Radian{(D - (1 + 2 * T1 + C1) * pow(D, 3) / 6 +
        (5 - 2 * C1 + 28 * T1 - 3 * pow(C1, 2) + 8 * EPrimeSquared + 24 * pow(T1, 2)) * pow(D, 5)
            / 120) /
        cos(phi1Rad)};
    longitude = Degree{lambda0 + Degree{longRad}};
}
SK42::SK42(WGS84 wgs_84) {
    altitude = wgs_84.altitude;
    latitude = Degree{wgs_84.latitude - dB(wgs_84.latitude, wgs_84.longitude, wgs_84.altitude, p) / 3600};
    longitude = Degree{wgs_84.longitude - dL(wgs_84.latitude, wgs_84.longitude, wgs_84.altitude, p) / 3600};
}

SK42::SK42(GaussKruger gk) {
    altitude = gk.height;

    int No = gk.y * pow(10, -6);
    double Bi = gk.x / 6367558.4968;
    double Bo = Bi + sin(Bi * 2) * (0.00252588685 - 0.0000149186 * pow(sin(Bi), 2) + 0.00000011904 * pow(sin(Bi), 4));
    double Zo = (gk.y - (10 * No + 5) * 100000) / (6378245 * cos(Bo));
    double Ba = Zo * Zo * (0.01672 - 0.0063 * pow(sin(Bo), 2) + 0.01188 * pow(sin(Bo), 4) - 0.00328 * pow(sin(Bo), 6));
    double Bb = Zo * Zo * (0.042858 - 0.025318 * pow(sin(Bo), 2) + 0.014346 * pow(sin(Bo), 4) - 0.001264 * pow(sin(Bo), 6) - Ba);
    double Bc = Zo * Zo * (0.10500614 - 0.04559916 * pow(sin(Bo), 2) + 0.00228901 * pow(sin(Bo), 4) - 0.00002987 * pow(sin(Bo), 6) - Bb);
    double dB = Zo * Zo * sin(Bo * 2) * (0.251684631 - 0.003369263 * pow(sin(Bo), 2) + 0.000011276 * pow(sin(Bo), 4) - Bc);
    latitude = Radian{Bo - dB};

    double La = Zo * Zo * (0.0038 + 0.0524 * pow(sin(Bo) , 2) + 0.0482 * pow(sin(Bo) , 4) + 0.0032 * pow(sin(Bo) , 6));
    double Lb = Zo * Zo * (0.01225 + 0.09477 * pow(sin(Bo) , 2) + 0.03282 * pow(sin(Bo) , 4) - 0.00034 * pow(sin(Bo) , 6) - La);
    double Lc = Zo * Zo * (0.0420025 + 0.1487407 * pow(sin(Bo) , 2) + 0.005942 * pow(sin(Bo) , 4) - 0.000015 * pow(sin(Bo) , 6) - Lb);
    double Ld = Zo * Zo * (0.16778975 + 0.16273586 * pow(sin(Bo) , 2) - 0.0005249 * pow(sin(Bo) , 4) - 0.00000846 * pow(sin(Bo) , 6) - Lc);
    double dL = Zo * (1 - 0.0033467108 * pow(sin(Bo) , 2) - 0.0000056002 * pow(sin(Bo) , 4) - 0.0000000187 * pow(sin(Bo) , 6) - Ld);
    longitude = Radian{Radian{Degree{6 * (No - 0.5)}} + dL};
}
GaussKruger::GaussKruger(SK42 sk_42) : height(sk_42.altitude) {
    double L = sk_42.longitude;
    Radian B = sk_42.latitude;
    int No = (6 + L) / 6;
    double Lo = Radian{Degree{L - (3 + 6 * (No - 1))}};
    double Xa = pow(Lo, 2) * (109500 - 574700 * pow(sin(B), 2) + 863700 * pow(sin(B), 4) - 398600 * pow(sin(B), 6));
    double Xb = pow(Lo, 2) * (278194 - 830174 * pow(sin(B), 2) + 572434 * pow(sin(B), 4) - 16010 * pow(sin(B), 6) + Xa);
    double Xc = pow(Lo, 2) * (672483.4 - 811219.9 * pow(sin(B), 2) + 5420 * pow(sin(B), 4) - 10.6 * pow(sin(B), 6) + Xb);
    double Xd = pow(Lo, 2) * (1594561.25 + 5336.535 * pow(sin(B), 2) + 26.79 * pow(sin(B), 4) + 0.149 * pow(sin(B), 6) + Xc);
    x = 6367558.4968 * B - sin(B * 2) * (16002.89 + 66.9607 * pow(sin(B), 2) + 0.3515 * pow(sin(B), 4) - Xd);

    double Ya = pow(Lo, 2) * (79690 - 866190 * pow(sin(B), 2) + 1730360 * pow(sin(B), 4) - 945460 * pow(sin(B), 6));
    double Yb = pow(Lo, 2) * (270806 - 1523417 * pow(sin(B), 2) + 1327645 * pow(sin(B), 4) - 21701 * pow(sin(B), 6) + Ya);
    double Yc = pow(Lo, 2) * (1070204.16 - 2136826.66 * pow(sin(B), 2) + 17.98 * pow(sin(B), 4) - 11.99 * pow(sin(B), 6) + Yb);
    y = (5 + 10 * No) * 100000 + Lo * cos(B) * (6378245 + 21346.1415 * pow(sin(B), 2) + 107.159 * pow(sin(B), 4) + 0.5977 * pow(sin(B), 6) + Yc);
}

PZ90::PZ90(Degree latitude, Degree longitude, double altitude)
    : latitude(latitude), longitude(longitude), altitude(altitude) {}
PZ90::PZ90(WGS84 wgs_84) {
    altitude = wgs_84.altitude;

    latitude = Degree{wgs_84.latitude - dB(wgs_84.latitude, wgs_84.longitude, wgs_84.altitude, p) / 3600};
    longitude = Degree{wgs_84.longitude - dL(wgs_84.latitude, wgs_84.longitude, wgs_84.altitude, p) / 3600};
}

char UTM::letter_designator(Degree latitude) {
    if ((84 >= latitude) && (latitude >= 72)) {
        return 'X';
    } else if ((72 > latitude) && (latitude >= 64)) {
        return 'W';
    } else if ((64 > latitude) && (latitude >= 56)) {
        return 'V';
    } else if ((56 > latitude) && (latitude >= 48)) {
        return 'U';
    } else if ((48 > latitude) && (latitude >= 40)) {
        return 'T';
    } else if ((40 > latitude) && (latitude >= 32)) {
        return 'S';
    } else if ((32 > latitude) && (latitude >= 24)) {
        return 'R';
    } else if ((24 > latitude) && (latitude >= 16)) {
        return 'Q';
    } else if ((16 > latitude) && (latitude >= 8)) {
        return 'P';
    } else if ((8 > latitude) && (latitude >= 0)) {
        return 'N';
    } else if ((0 > latitude) && (latitude >= -8)) {
        return 'M';
    } else if ((-8 > latitude) && (latitude >= -16)) {
        return 'L';
    } else if ((-16 > latitude) && (latitude >= -24)) {
        return 'K';
    } else if ((-24 > latitude) && (latitude >= -32)) {
        return 'J';
    } else if ((-32 > latitude) && (latitude >= -40)) {
        return 'H';
    } else if ((-40 > latitude) && (latitude >= -48)) {
        return 'G';
    } else if ((-48 > latitude) && (latitude >= -56)) {
        return 'F';
    } else if ((-56 > latitude) && (latitude >= -64)) {
        return 'E';
    } else if ((-64 > latitude) && (latitude >= -72)) {
        return 'D';
    } else if ((-72 > latitude) && (latitude >= -80)) {
        return 'C';
    }
    // 'Z' is an error flag, the latitudeitude is outside the UTM limits
    return 'Z';
}
UTM::UTM(Degree E, Degree N, double altitude, std::string  zone)
    : E(E), N(N), altitude(altitude), zone(std::move(zone)) {}
UTM::UTM(WGS84 wgs_84) {
    altitude = wgs_84.altitude;

    Radian latRad = wgs_84.latitude;
    Radian longRad = wgs_84.longitude;

    int zoneNumber = static_cast<int>((wgs_84.longitude + 180) / 6) + 1;

    if (wgs_84.latitude >= 56.0 && wgs_84.latitude < 64.0 && wgs_84.longitude >= 3.0 && wgs_84.longitude < 12.0) {
        zoneNumber = 32;
    }

    // Special zones for Svalbard
    if (wgs_84.latitude >= 72.0 && wgs_84.latitude < 84.0) {
        if (wgs_84.longitude >= 0.0 && wgs_84.longitude < 9.0) {
            zoneNumber = 31;
        } else if (wgs_84.longitude >= 9.0 && wgs_84.longitude < 21.0) {
            zoneNumber = 33;
        } else if (wgs_84.longitude >= 21.0 && wgs_84.longitude < 33.0) {
            zoneNumber = 35;
        } else if (wgs_84.longitude >= 33.0 && wgs_84.longitude < 42.0) {
            zoneNumber = 37;
        }
    }
    // +3 puts origin in middle of zone
    Degree lambda0{(zoneNumber - 1) * 6 - 177};
    Radian lambda0Rad = lambda0;

    // Compute the UTM Zone from the latitude and longitude
    zone = std::to_string(zoneNumber) + letter_designator(wgs_84.latitude);

    constexpr double EPrimeSquared = WGS84::_e2 / (1 - WGS84::_e2);

    double N_ = WGS84::_a / sqrt(1 - WGS84::_e2 * pow(sin(latRad), 2));
    double T = tan(latRad) * tan(latRad);
    double C = EPrimeSquared * pow(cos(latRad), 2);
    double A = cos(latRad) * (longRad - lambda0Rad);

    double M = WGS84::_a * latRad * ((1 - WGS84::_e2 / 4 - (3. / 64) * pow(WGS84::_e2, 2) - (5. / 256) * pow(WGS84::_e2, 3)) -
            sin(2 * latRad) * ((3. / 8) * WGS84::_e2 + (3. / 32) * pow(WGS84::_e2, 2) + (45. / 1024) * pow(WGS84::_e2, 3)) +
            sin(4 * latRad) * ((15. / 256) * pow(WGS84::_e2, 2) + (45. / 1024) * pow(WGS84::_e2, 3)) -
            sin(6 * latRad) * ((35. / 3072) * pow(WGS84::_e2, 3)));

    E = k0 * N_ * (A + (1 - T + C) * pow(A, 3) / 6 +
        (5 - 18 * T + pow(T, 2) + 72 * C - 58 * EPrimeSquared) * pow(A, 5) / 120) + E0;

    N = k0 * (M + N_ * tan(latRad) *
        (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * pow(A, 4) / 24 +
            (61 - 58 * T + T * T + 600 * C - 330 * EPrimeSquared) * pow(A, 5) / 720));

    if (wgs_84.latitude < 0) {
        // 10000000 meter offset for southern hemisphere
        N += N0;
    }
}
