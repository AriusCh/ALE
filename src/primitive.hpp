#ifndef ALE_SOLVER_SRC_PRIMITIVE_HPP_
#define ALE_SOLVER_SRC_PRIMITIVE_HPP_

#include <array>
#include <memory>

#include "boundary.hpp"

enum class PrimitiveType {
  ePoint,
  eLine,
  eEdge,
  ePolygon,
  eRectangle,
  eCircle
};

class Point {
 public:
  Point(double x_, double y_);
  Point(const Point &rhs) = default;
  Point(Point &&rhs) = default;

  Point &operator=(const Point &rhs) = default;
  Point &operator=(Point &&rhs) = default;

  ~Point() = default;

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
  Line(Point first_, Point second_);

  Point getFirstPoint() const;
  Point getSecondPoint() const;

  void setFirstPoint(Point point);
  void setSecondPoint(Point point);

  double getLength() const;

 private:
  Point first;
  Point second;

  PrimitiveType type = PrimitiveType::eLine;
};

class Edge {
 public:
  Edge(const std::vector<Line> &lines_);
  Edge(const Edge &rhs) = default;
  Edge(Edge &&rhs) = default;

  Edge &operator=(const Edge &rhs) = default;
  Edge &operator=(Edge &&rhs) = default;

 public:
  const std::vector<Line> &getLines() const;
  std::vector<Line> &getLines();

  bool isContinuous() const;

 private:
  std::vector<Line> lines;

  PrimitiveType type = PrimitiveType::eEdge;
};

class Polygon {
 public:
  Polygon(std::shared_ptr<Edge> leftEdge, std::shared_ptr<Edge> topEdge,
          std::shared_ptr<Edge> rightEdge, std::shared_ptr<Edge> bottomEdge);

 public:
  bool isEnclosed() const;

  virtual bool isInside(double x_, double y_) const;

 private:
  bool isLineIntersecting(double x_, double y_, const Line &line) const;

 protected:
  std::shared_ptr<Edge> leftEdge;
  std::shared_ptr<Edge> topEdge;
  std::shared_ptr<Edge> rightEdge;
  std::shared_ptr<Edge> bottomEdge;
};

class Rectangle : public Polygon {
 public:
  Rectangle(double x0, double y0, double width, double height);

 public:
  virtual bool isInside(double x_, double y_) const override;
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
