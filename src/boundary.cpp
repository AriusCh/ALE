#include "boundary.hpp"
Boundary::Boundary(BoundaryType type_) : type(type_) {}

ExternalBoundary::ExternalBoundary(BoundaryType type_,
                                   ExternalBoundarySide side_)
    : Boundary(type_), side(side_) {}
