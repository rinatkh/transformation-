//
// Created by kira on 07.09.2020.
//

#include "radian_degree.h"

#include <cmath>

Degree::Degree(Radian radian) {
    _degree = radian * 180 / M_PI;
}
Degree::operator double() const {
    return _degree;
}

Radian::Radian(Degree degree) {
    _radian = degree * M_PI / 180;
}
Radian::operator double() const {
    return _radian;
}
