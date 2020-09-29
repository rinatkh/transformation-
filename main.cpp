#include "transformations.h"

#include <iostream>

int main() {
    enum COMMAND {
        FromWGS84ToGaussKruger = 1,
        FromGaussKrugerToWGS84 = 2,
        FromPZ90ToGaussKruger = 3,
        FromGaussKrugerToPZ90 = 4,
        FromWGS84ToUTM = 5,
        FromUTMToWGS84 = 6,
        FromPZ90ToUTM = 7,
        FromUTMToPZ90 = 8
    };
    std::cout << "enter command:\n"
                 "1: FromWGS84ToGaussKruger\n"
                 "2: FromGaussKrugerToWGS84\n"
                 "3: FromPZ90ToGaussKruger\n"
                 "4: FromGaussKrugerToPZ90\n"
                 "5: FromWGS84ToUTM\n"
                 "6: FromUTMToWGS84\n"
                 "7: FromPZ90ToUTM\n"
                 "8: FromUTMToPZ90"
              << std::endl;
    int command = 0;
    std::cout.precision(9);
    while (std::cin >> command) {
        switch (command) {
            case FromWGS84ToGaussKruger: {
                std::cout << "latitude: ";
                double latitude{};
                std::cin >> latitude;

                std::cout << "longitude: ";
                double longitude{};
                std::cin >> longitude;

                std::cout << "altitude: ";
                double altitude{};
                std::cin >> altitude;

                WGS84 wgs_84{Degree{latitude}, Degree{longitude}, altitude};
                GaussKruger gk{SK42{wgs_84}};

                std::cout << "x: " << gk.x << " y: " << gk.y << std::endl;
                break;
            }
            case FromGaussKrugerToWGS84: {
                GaussKruger gk{};
                std::cout << "x: ";
                std::cin >> gk.x;
                std::cout << "y: ";
                std::cin >> gk.y;
                std::cout << "height: ";
                std::cin >> gk.height;

                WGS84 wgs_84{SK42{gk}};

                std::cout << "latitude: " << wgs_84.latitude << " longitude: " << wgs_84.longitude << std::endl;
                break;
            }
            case FromPZ90ToGaussKruger: {
                std::cout << "latitude: ";
                double latitude{};
                std::cin >> latitude;

                std::cout << "longitude: ";
                double longitude{};
                std::cin >> longitude;

                std::cout << "altitude: ";
                double altitude{};
                std::cin >> altitude;

                PZ90 pz_90{Degree{latitude}, Degree{longitude}, altitude};
                GaussKruger gk{SK42{WGS84{pz_90}}};

                std::cout << "x: " << gk.x << " y: " << gk.y << std::endl;
                break;
            }
            case FromGaussKrugerToPZ90: {
                GaussKruger gk{};
                std::cout << "x: ";
                std::cin >> gk.x;
                std::cout << "y: ";
                std::cin >> gk.y;
                std::cout << "height: ";
                std::cin >> gk.height;

                PZ90 pz_90{WGS84{SK42{gk}}};

                std::cout << "latitude: " << pz_90.latitude << " longitude: " << pz_90.longitude << std::endl;
                break;
            }
            case FromWGS84ToUTM: {
                std::cout << "latitude: ";
                double latitude{};
                std::cin >> latitude;

                std::cout << "longitude: ";
                double longitude{};
                std::cin >> longitude;

                std::cout << "altitude: ";
                double altitude{};
                std::cin >> altitude;

                WGS84 wgs_84{Degree{latitude}, Degree{longitude}, altitude};
                UTM utm{wgs_84};

                std::cout << "Easting: " << utm.E << " Northing: " << utm.N << " Zone: " << utm.zone
                          << std::endl;
                break;
            }
            case FromUTMToWGS84: {
                std::cout << "Easting: ";
                double E{};
                std::cin >> E;

                std::cout << "Northing: ";
                double N{};
                std::cin >> N;

                std::cout << "altitude: ";
                double altitude{};
                std::cin >> altitude;

                std::cout << "Zone: ";
                std::string zone;
                std::cin >> zone;

                UTM utm{Degree{E}, Degree{N}, altitude, zone};
                WGS84 wgs_84{utm};

                std::cout << "latitude: " << wgs_84.latitude << " longitude: " << wgs_84.longitude << std::endl;
                break;
            }
            case FromPZ90ToUTM: {
                std::cout << "latitude: ";
                double latitude{};
                std::cin >> latitude;

                std::cout << "longitude: ";
                double longitude{};
                std::cin >> longitude;

                std::cout << "altitude: ";
                double altitude{};
                std::cin >> altitude;

                PZ90 pz_90{Degree{latitude}, Degree{longitude}, altitude};
                UTM utm{WGS84{pz_90}};

                std::cout << "Easting: " << utm.E << " Northing: " << utm.N << " Zone: " << utm.zone
                          << std::endl;
                break;
            }
            case FromUTMToPZ90: {
                std::cout << "Easting: ";
                double E{};
                std::cin >> E;

                std::cout << "Northing: ";
                double N{};
                std::cin >> N;

                std::cout << "altitude: ";
                double altitude{};
                std::cin >> altitude;

                std::cout << "Zone: ";
                std::string zone;
                std::cin >> zone;

                UTM utm{Degree{E}, Degree{N}, altitude, zone};
                PZ90 pz_90{WGS84{utm}};

                std::cout << "latitude: " << pz_90.latitude << " longitude: " << pz_90.longitude << std::endl;
                break;
            }
            default:
                std::cout << "wrong command" << std::endl;
        }
    }
    return 0;
}
