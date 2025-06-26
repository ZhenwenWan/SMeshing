#include "geometricModel.hpp"
#include <catch.hpp>

#include <iostream>

TEST_CASE("geometricModel", "[geometricModel]")
{
    using namespace shonCloud;

    SECTION("boundingBox")
    {
        // auto thisGeoModel = std::make_shared<geometricModel>();
        // const std::string STEPFileName("testData/geometries/box.step");
        // thisGeoModel->readSTEP(STEPFileName);
        // REQUIRE(thisGeoModel->boundingBox_.first[0] == Approx(-0.0000001));
        // REQUIRE(thisGeoModel->boundingBox_.first[1] == Approx(-0.0000001));
        // REQUIRE(thisGeoModel->boundingBox_.first[2] == Approx(-0.0000001));
        // REQUIRE(thisGeoModel->boundingBox_.second[0] == Approx(0.2));
        // REQUIRE(thisGeoModel->boundingBox_.second[1] == Approx(0.2));
        // REQUIRE(thisGeoModel->boundingBox_.second[2] == Approx(0.2));
        // REQUIRE(thisGeoModel->solids_.size() == 1);
        // REQUIRE(thisGeoModel->edges_.size() == 12);
        // REQUIRE(thisGeoModel->wires_.size() == 6);
        // REQUIRE(thisGeoModel->vertices_.size() == 8);
        // REQUIRE(thisGeoModel->faces_.size() == 6);
        // for (const auto& oneSolid : thisGeoModel->solids_)
        // {
        //     std::cout << " solid tag " << oneSolid.second->tag_ << std::endl;
        //     for (const auto& oneFace : oneSolid.second->faces_)
        //     {
        //         std::cout << " face tag: " << oneFace.second->tag_ <<
        //         std::endl; for (const auto& oneWire : oneFace.second->wires_)
        //         {
        //             std::cout << "  wire tag: " << oneWire.second->tag_
        //                       << std::endl;
        //             for (const auto& oneEdge : oneWire.second->edges_)
        //             {
        //                 std::cout << "   edge tag: " << oneEdge.second->tag_
        //                           << std::endl;

        //                 for (const auto& oneVertex :
        //                 oneEdge.second->vertices_)
        //                 {
        //                     std::cout
        //                         << "    vertex tag: " <<
        //                         oneVertex.second->tag_
        //                         << std::endl;
        //                 }
        //             }
        //         }
        //     }
        // }
    }
    SECTION("edgeParameter")
    {
        auto thisGeoModel = std::make_shared<geometricModel>();
        const std::string STEPFileName("testData/geometries/cutBox.step");
        thisGeoModel->readSTEP(STEPFileName);
        // for (const auto& oneEdge : thisGeoModel->edges_)
        // {
        //     std::cout << oneEdge.second->uBounds_.first << " "
        //               << oneEdge.second->uBounds_.second << std::endl;
        // }
    }
}
