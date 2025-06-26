#include "tetMesh.hpp"
#include <catch.hpp>

#include "shonModelerEnums.hpp"

#include "boost/filesystem.hpp"
#include "shonModelerEnums.hpp"
#include <iostream>

TEST_CASE("tetMesh", "[tetMesh]")
{
    using namespace shonCloud;

    // SECTION("fromSurface")
    //{
    //    const std::string stlFileName("testData/geometries/cutBox.stl");
    //    auto thisTetMesh = std::make_shared<tetMesh>();
    //    thisTetMesh->generate3DMesh(stlFileName);
    //}

    // SECTION("basicTest")
    // {
    //     std::default_random_engine eng(std::random_device{}());
    //     std::uniform_real_distribution<double> dist_w(0, 1);
    //     std::uniform_real_distribution<double> dist_h(0, 1);

    //     const UIN numberOfPoints = 40;

    //     auto thisTetMesh = std::make_shared<tetMesh>();

    //     for (UIN i = 0; i < numberOfPoints; ++i)
    //     {
    //         vec3d oneVec(dist_w(eng), dist_h(eng), 0);
    //         auto oneNode = std::make_shared<node>(oneVec);
    //         thisTetMesh->nodes_.push_back(oneNode);
    //     }

    //     thisTetMesh->triangulateBowyerWatson();
    //     const std::string fileName = "/home/martin/Documents/test/mesh.vtu";
    //     thisTetMesh->writeVTK(fileName);
    //     REQUIRE(2 == 2);
    // }

    // SECTION("fromSurface")
    // {
    //     const std::string stlFileName("testData/geometries/cutBox.stl");
    //     auto thisTetMesh = std::make_shared<tetMesh>();
    //     thisTetMesh->generate2DMesh(stlFileName);
    // }

    // SECTION("triangulate3d")
    // {
    //     const std::string stlFileName("testData/geometries/cutBox.stl");
    //     auto thisTetMesh = std::make_shared<tetMesh>();
    //     thisTetMesh->generate3DMesh(stlFileName);

    //     const std::string fileName = "/home/martin/Documents/test/mesh.vtu";
    //     thisTetMesh->writeVTKTetra(fileName);
    // }

    SECTION("seed2DPoints_CutBox")
    {
        const std::string cwd = boost::filesystem::current_path().string();
        auto thisTetMesh = std::make_shared<tetMesh>();
        const std::string STEPFileName("testData/geometries/cutBox.step");
        thisTetMesh->readFromStep(STEPFileName);
        // for (const auto& oneSolid :
        // thisTetMesh->geometricModel_->solids_)
        thisTetMesh->seedVerticesOnEdges(0.01);
        for (const auto& oneSolid : thisTetMesh->geometricModel_->faces_)
        {
            // if (oneSolid.second->tag_ > 0) continue;
            vectorVec3d nodes;
            // oneSolid.second->tetraBowyerWatson();
            oneSolid.second->mesh(0.01);
            const std::string resultPath =
                cwd +
                "/testData/results/"
                "1DCutBox-Face-" +
                std::to_string(oneSolid.second->tag_) + ".vtu";
            oneSolid.second->writeVerticesToVTK(resultPath);
            const std::string resultPath1 =
                cwd +
                "/testData/results/"
                "CutBox-Solid-" +
                std::to_string(oneSolid.second->tag_) + ".vtu";
            oneSolid.second->writeVTK(resultPath1);
        }
    }
}
