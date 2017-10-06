///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012-2017 DreamWorks Animation LLC
//
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
//
// Redistributions of source code must retain the above copyright
// and license notice and the following restrictions and disclaimer.
//
// *     Neither the name of DreamWorks Animation nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// IN NO EVENT SHALL THE COPYRIGHT HOLDERS' AND CONTRIBUTORS' AGGREGATE
// LIABILITY FOR ALL CLAIMS REGARDLESS OF THEIR BASIS EXCEED US$250.00.
//
///////////////////////////////////////////////////////////////////////////

#include <random>

#include <cppunit/extensions/HelperMacros.h>

#include <openvdb/openvdb.h>

#include <openvdb/points/PointScatter.h>
#include <openvdb/points/PointCount.h>
#include <openvdb/points/PointDataGrid.h>

#include <openvdb/math/Math.h>
#include <openvdb/math/Coord.h>

using namespace openvdb;
using namespace openvdb::points;

class TestPointScatter: public CppUnit::TestCase
{
public:

    virtual void setUp() { openvdb::initialize(); }
    virtual void tearDown() { openvdb::uninitialize(); }

    CPPUNIT_TEST_SUITE(TestPointScatter);
    CPPUNIT_TEST(testScatter);
    CPPUNIT_TEST(testVoxelScatter);
    CPPUNIT_TEST(testWeightedVoxelScatter);
    CPPUNIT_TEST_SUITE_END();

    void testScatter();
    void testVoxelScatter();
    void testWeightedVoxelScatter();

}; // class TestPointScatter

void
TestPointScatter::testScatter(){}

void
TestPointScatter::testVoxelScatter()
{
    const Index32 pointsPerVoxel = 8;
    const math::CoordBBox boxBounds(math::Coord(0), math::Coord(1)); // 8 voxels

    // Test the free function for all default grid types

    {
        BoolGrid grid;
        grid.sparseFill(boxBounds, false, /*active*/true);
        auto points = points::pointVoxelScatter
            <BoolGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        DoubleGrid grid;
        grid.sparseFill(boxBounds, 0.0, /*active*/true);
        auto points = points::pointVoxelScatter
            <DoubleGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        FloatGrid grid;
        grid.sparseFill(boxBounds, 0.0f, /*active*/true);
        auto points = points::pointVoxelScatter
            <FloatGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        Int32Grid grid;
        grid.sparseFill(boxBounds, 0, /*active*/true);
        auto points = points::pointVoxelScatter
            <Int32Grid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        Int64Grid grid;
        grid.sparseFill(boxBounds, 0, /*active*/true);
        auto points = points::pointVoxelScatter
            <Int64Grid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        MaskGrid grid;
        grid.sparseFill(boxBounds, /*maskBuffer*/true);
        auto points = points::pointVoxelScatter
            <MaskGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        StringGrid grid;
        grid.sparseFill(boxBounds, "", /*active*/true);
        auto points = points::pointVoxelScatter
            <StringGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        Vec3DGrid grid;
        grid.sparseFill(boxBounds, Vec3d(), /*active*/true);
        auto points = points::pointVoxelScatter
            <Vec3DGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        Vec3IGrid grid;
        grid.sparseFill(boxBounds, Vec3i(), /*active*/true);
        auto points = points::pointVoxelScatter
            <Vec3IGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        Vec3SGrid grid;
        grid.sparseFill(boxBounds, Vec3f(), /*active*/true);
        auto points = points::pointVoxelScatter
            <Vec3SGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }
    {
        PointDataGrid grid;
        grid.sparseFill(boxBounds, 0, /*active*/true);
        auto points = points::pointVoxelScatter
            <PointDataGrid, std::mt19937, PointDataGrid>(grid, 0, pointsPerVoxel);
        CPPUNIT_ASSERT_EQUAL(Index32(1), points->tree().leafCount());
        CPPUNIT_ASSERT_EQUAL(Index64(8), points->activeVoxelCount());
        CPPUNIT_ASSERT_EQUAL(Index64(pointsPerVoxel * 8), pointCount(points->tree()));
    }

    std::mt19937 engine19937;
    UniformVoxelPointScatter<std::mt19937> scatter19937(engine19937, pointsPerVoxel, /*seed=*/0);

    // Test a grid containing tiles scatters correctly

    BoolGrid grid;
    grid.tree().addTile(/*level*/1, math::Coord(0), /*value*/true, /*active*/true);
    grid.tree().setValueOn(math::Coord(8,0,0)); // add another leaf

    const Index32 NUM_VALUES = BoolGrid::TreeType::LeafNodeType::NUM_VALUES;

    CPPUNIT_ASSERT_EQUAL(Index32(1), grid.tree().leafCount());
    CPPUNIT_ASSERT_EQUAL(Index64(NUM_VALUES + 1), grid.activeVoxelCount());

    scatter19937(grid);
    auto points = scatter19937.points();

    const Index64 expectedCount = Index64(pointsPerVoxel * (NUM_VALUES + 1));

#ifndef OPENVDB_2_ABI_COMPATIBLE
    CPPUNIT_ASSERT_EQUAL(Index64(0), points->tree().activeTileCount());
#endif
    CPPUNIT_ASSERT_EQUAL(Index32(2), points->tree().leafCount());
    CPPUNIT_ASSERT_EQUAL(Index64(NUM_VALUES + 1), points->activeVoxelCount());
    CPPUNIT_ASSERT_EQUAL(expectedCount, pointCount(points->tree()));

    // Explicitly check P attribute

    const auto* attributeSet = &(points->tree().cbeginLeaf()->attributeSet());
    CPPUNIT_ASSERT_EQUAL(size_t(1), attributeSet->size());
    const auto* array = attributeSet->getConst(0);
    CPPUNIT_ASSERT(array);

    using PositionArrayT = TypedAttributeArray<Vec3f, NullCodec>;
    CPPUNIT_ASSERT(array->isType<PositionArrayT>());

    size_t size = array->size();
    CPPUNIT_ASSERT_EQUAL(size_t(pointsPerVoxel * NUM_VALUES), size);

    AttributeHandle<Vec3f, NullCodec>::Ptr pHandle =
        AttributeHandle<Vec3f, NullCodec>::create(*array);
    for (size_t i = 0; i < size; ++i) {
        const Vec3f P = pHandle->get(i);
        CPPUNIT_ASSERT(P[0] >=-0.5f);
        CPPUNIT_ASSERT(P[0] <= 0.5f);
        CPPUNIT_ASSERT(P[1] >=-0.5f);
        CPPUNIT_ASSERT(P[1] <= 0.5f);
        CPPUNIT_ASSERT(P[2] >=-0.5f);
        CPPUNIT_ASSERT(P[2] <= 0.5f);
    }

    // Test the rng seed

    const Vec3f firstPosition = pHandle->get(0);
    scatter19937.setSeed(1);
    scatter19937(grid);
    points = scatter19937.points();

    attributeSet = &(points->tree().cbeginLeaf()->attributeSet());
    CPPUNIT_ASSERT_EQUAL(size_t(1), attributeSet->size());

    array = attributeSet->getConst(0);
    CPPUNIT_ASSERT(array);
    CPPUNIT_ASSERT(array->isType<PositionArrayT>());

    size = array->size();
    CPPUNIT_ASSERT_EQUAL(size_t(pointsPerVoxel * NUM_VALUES), size);
    pHandle = AttributeHandle<Vec3f, NullCodec>::create(*array);

    const Vec3f secondPosition = pHandle->get(0);
    CPPUNIT_ASSERT(firstPosition[0] != secondPosition[0]);
    CPPUNIT_ASSERT(firstPosition[1] != secondPosition[1]);
    CPPUNIT_ASSERT(firstPosition[2] != secondPosition[2]);

    // Test spread

    scatter19937.setSpread(0.2f);
    scatter19937(grid);
    points = scatter19937.points();

    attributeSet = &(points->tree().cbeginLeaf()->attributeSet());
    CPPUNIT_ASSERT_EQUAL(size_t(1), attributeSet->size());
    array = attributeSet->getConst(0);
    CPPUNIT_ASSERT(array);
    CPPUNIT_ASSERT(array->isType<PositionArrayT>());

    size = array->size();
    CPPUNIT_ASSERT_EQUAL(size_t(pointsPerVoxel * NUM_VALUES), size);

    pHandle = AttributeHandle<Vec3f, NullCodec>::create(*array);
    for (size_t i = 0; i < size; ++i) {
        const Vec3f P = pHandle->get(i);
        CPPUNIT_ASSERT(P[0] >=-0.2f);
        CPPUNIT_ASSERT(P[0] <= 0.2f);
        CPPUNIT_ASSERT(P[1] >=-0.2f);
        CPPUNIT_ASSERT(P[1] <= 0.2f);
        CPPUNIT_ASSERT(P[2] >=-0.2f);
        CPPUNIT_ASSERT(P[2] <= 0.2f);
    }

    // Test mt11213b

    using mt11213b = std::mersenne_twister_engine<uint32_t, 32, 351, 175, 19,
        0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253>;

    mt11213b engine11213b;
    UniformVoxelPointScatter<mt11213b> scatter11213b(engine11213b, pointsPerVoxel, /*seed=*/0);

    scatter11213b(grid);
    points = scatter11213b.points();

    CPPUNIT_ASSERT_EQUAL(Index32(2), points->tree().leafCount());
    CPPUNIT_ASSERT_EQUAL(Index64(NUM_VALUES + 1), points->activeVoxelCount());
    CPPUNIT_ASSERT_EQUAL(expectedCount, pointCount(points->tree()));
}

void
TestPointScatter::testWeightedVoxelScatter(){}


CPPUNIT_TEST_SUITE_REGISTRATION(TestPointScatter);

// Copyright (c) 2012-2017 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
