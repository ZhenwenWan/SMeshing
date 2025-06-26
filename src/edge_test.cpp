#include "edge.hpp"
#include <catch.hpp>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <gp_Pnt.hxx>

#include <iostream>

TEST_CASE("edge", "[edge]")
{
    using namespace shonCloud;
    SECTION("parameterTransformation")
    {
        vec3d p1(0, 0, 0);
        vec3d p2(0.2, 0, 0);
        gp_Pnt gp1(p1[0], p1[1], p1[2]);
        gp_Pnt gp2(p2[0], p2[1], p2[2]);
        TopoDS_Edge occEdge = BRepBuilderAPI_MakeEdge(gp1, gp2).Edge();
        int edgeTag = 0;
        auto oneEdge = std::make_shared<edge>(occEdge, edgeTag);
        REQUIRE(oneEdge->uBounds_.first == 0.0);
        REQUIRE(oneEdge->uBounds_.second == 0.2);
        const vec3d onePoint(0.1, 0, 0);
        const double u1 = oneEdge->parameterFromPoint(onePoint);
        REQUIRE(u1 == 0.1);
        REQUIRE(oneEdge->pointFromParameter(u1) == onePoint);
    }
}
