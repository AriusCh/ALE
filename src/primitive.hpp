#ifndef ALE_SOLVER_SRC_PRIMITIVE_HPP_
#define ALE_SOLVER_SRC_PRIMITIVE_HPP_

#include <array>
#include <memory>

#include "boundary.hpp"

enum class PrimitiveType { ePoint, eLine, eEdge, eRectangle, eCircle };

class Point {
 public:
  Point(double x_, double y_);

  bool operator==(const Point &) const = default;
 public:
  double getX() const;
  double getY() const;

  void setX(double x_);
  void setY(double y_);

 private:
  double x;
  double y;

  PrimitiveType type = PrimitiveType::ePoint;
};

class Line {
 public:
  Line(std::shared_ptr<Point> first_, std::shared_ptr<Point> second_);

  std::shared_ptr<Point> getFirstPoint() const;
  std::shared_ptr<Point> getSecondPoint() const;

  void setFirstPoint(std::shared_ptr<Point> point);
  void setSecondPoint(std::shared_ptr<Point> point);

  double getLength() const;

 private:
  std::shared_ptr<Point> first;
  std::shared_ptr<Point> second;

  PrimitiveType type = PrimitiveType::eLine;
};

class Edge {
 public:
  Edge(const std::vector<Line> &lines);

 public:
  const std::vector<Line> &getLines() const;
  std::vector<Line> &getLines();

  bool isContinuous() const;

 private:
  std::vector<Line> lines;

  PrimitiveType type = PrimitiveType::eEdge;
};

class Polygon {
  Polygon();

};

// class Primitive2D {
//  public:
//   virtual bool isInside(double x, double y) const;
// };
//
// class PrimitiveRectangle : public Primitive2D {
//  public:
//   PrimitiveRectangle(double x0, double y0, double width, double height,
//                      double phi = 0.);
//
//  public:
//   virtual bool isInside(double x_, double y_) const override;
//
//   void getVertices(std::array<double, 4> &x, std::array<double, 4> &y) const;
//
//   void setBoundaries(const std::vector<BoundaryType> &bnds);
//   std::vector<BoundaryType> getBoundaries() const;
//
//  private:
//   std::array<double, 4> x;
//   std::array<double, 4> y;
//
//   std::array<BoundaryType, 4> boundaries;
// };
//
// class PrimitiveCircle : public Primitive2D {
//  public:
//   PrimitiveCircle(double x0_, double y0_, double r);
//
//  public:
//   virtual bool isInside(double x_, double y_) const override;
//
//  private:
//   double x0;
//   double y0;
//   double r;
//
//   BoundaryType boundary;
// };
//
#endif
