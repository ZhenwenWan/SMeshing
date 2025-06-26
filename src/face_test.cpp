#include "face.hpp"
#include <catch.hpp>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRep_Builder.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>

#include <iostream>

TEST_CASE("occFace", "[occFace]")
{
    using namespace shonCloud;
    SECTION("boundingBox")
    {
        BRep_Builder aBuilder;
        gp_Pln planeXY;
        TopoDS_Face aFace = BRepBuilderAPI_MakeFace(planeXY);
        gp_Ax2 Ax2(gp_Pnt(2.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0),
                   gp_Dir(1.0, 0.0, 0.0));
        TopoDS_Wire wireOut =
            BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(gp_Circ(Ax2, 1.5)));
        aBuilder.Add(aFace, wireOut);
        const int faceTag = 0;
        auto oneFace = std::make_shared<face>(aFace, faceTag);
        REQUIRE(oneFace->boundingBox_.first[0] == Approx(0.4999999));
        REQUIRE(oneFace->boundingBox_.first[1] == Approx(-1.5000001));
        REQUIRE(oneFace->boundingBox_.first[2] == Approx(-0.0000001));
        REQUIRE(oneFace->boundingBox_.second[0] == Approx(3.5000001));
        REQUIRE(oneFace->boundingBox_.second[1] == Approx(1.5));
        REQUIRE(oneFace->boundingBox_.second[2] == Approx(0.0000001));
    }
}
