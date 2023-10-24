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

/* Point is a dot with coordinates (x, y) */
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

/* Line is a directional vector with start at first Point and end at the second
 */
class Line {
 public:
  Line(Point first_, Point second_);
  Line(const Line &rhs) = default;
  Line(Line &&rhs) = default;

  Line &operator=(const Line &rhs) = default;

 public:
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

/* Edge is an array of continuous lines which will correspond to a particular
 * side of a quadrilateral grid */
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

/* Polygon is an enclosed space which contains edges */
class Polygon {
 public:
  Polygon(std::shared_ptr<Edge> leftEdge, std::shared_ptr<Edge> topEdge,
          std::shared_ptr<Edge> rightEdge, std::shared_ptr<Edge> bottomEdge);

 public:
  virtual void generateMesh(std::vector<std::vector<double>> &x,
                            std::vector<std::vector<double>> &y) const;

  bool isEnclosed() const;
  virtual bool isInside(double x_, double y_) const;

  std::shared_ptr<Edge> getLeftEdge() const;
  std::shared_ptr<Edge> getTopEdge() const;
  std::shared_ptr<Edge> getRightEdge() const;
  std::shared_ptr<Edge> getBottomEdge() const;
  PrimitiveType getType() const;

 private:
  bool isLineIntersecting(double x_, double y_, const Line &line) const;

 protected:
  std::shared_ptr<Edge> leftEdge;
  std::shared_ptr<Edge> topEdge;
  std::shared_ptr<Edge> rightEdge;
  std::shared_ptr<Edge> bottomEdge;

 private:
  PrimitiveType type;
};

/* rectangle is a polygon with four straight perpendicular edges */
class Rectangle : public Polygon {
 public:
  Rectangle(double x0, double y0, double width, double height);
  Rectangle(double x0, double y0, double width, double height, double phi);

 public:
  virtual void generateMesh(std::vector<std::vector<double>> &x,
                            std::vector<std::vector<double>> &y) const override;

  virtual bool isInside(double x_, double y_) const override;
};

/* Region2D is an array of polygons to support arbitrary shape */
class Region2D {
 public:
  Region2D(std::unique_ptr<Polygon> &&polygon);
  Region2D(std::vector<std::unique_ptr<Polygon>> &&polygons_);

 public:
  bool isInside(double x_, double y_) const;

  const std::vector<std::unique_ptr<Polygon>> &getPolygons() const;
  std::vector<std::unique_ptr<Polygon>> &getPolygons();

 private:
  std::vector<std::unique_ptr<Polygon>> polygons;
};

#endif
