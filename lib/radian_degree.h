//
// Created by kira on 07.09.2020.
//

#ifndef MAIN_LIB_RADIAN_DEGREE_H_
#define MAIN_LIB_RADIAN_DEGREE_H_

class Degree;
class Radian;

class Degree {
 public:
    Degree() = default;
    explicit Degree(int degree) : _degree(degree) {}
    explicit Degree(double degree) : _degree(degree) {}
    Degree(Radian radian);
    operator double() const;

 private:
    double _degree;
};

class Radian {
 public:
    Radian() = default;
    explicit Radian(double radian) : _radian(radian) {}
    Radian(Degree degree);
    operator double() const;

 private:
    double _radian;
};


#endif  // MAIN_LIB_RADIAN_DEGREE_H_
