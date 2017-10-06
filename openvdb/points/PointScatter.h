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
//
/// @author Nick Avramoussis
///
/// @file PointScatter.h
///
/// @brief Various point scattering methods for generating VDB Points.
///
///  Note that this architecture is heavily based off the similar
///  functionality for more generic position scattering in
///  tools/PointScatter.h, with an additional set of free functions for
///  convenience.
///
#ifndef OPENVDB_POINTS_POINT_SCATTER_HAS_BEEN_INCLUDED
#define OPENVDB_POINTS_POINT_SCATTER_HAS_BEEN_INCLUDED

#include <algorithm>

#include <openvdb/openvdb.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/util/NullInterrupter.h>

#include "PointDataGrid.h"
#include "AttributeArray.h"
#include "PointCount.h"

#include <tbb/parallel_sort.h>
#include <tbb/parallel_for.h>

namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace points {

/// @brief Uniformly scatter a total amount of points in active regions
///
/// @param grid   A source grid. The resulting PointDataGrid will copy this grids
///               transform and voxelized topology.
/// @param seed   A seed for the RandGenT
/// @param count  The total number of points to scatter
/// @note returns the scattered PointDataGrid
template<
    typename GridT,
    typename RandGenT,
    typename PointDataGridT>
inline typename PointDataGridT::Ptr
pointScatter(const GridT& grid,
             const unsigned int seed,
             const Index64 count);

/// @brief Uniformly scatter a fixed number of points per active voxel. If the pointsPerVoxel
///        value provided is a fractional value, each voxel calculates a delta value of
///        how likely it is to contain an extra point.
///
/// @param grid            A source grid. The resulting PointDataGrid will copy this grids
///                        transform and voxelized topology.
/// @param seed            A seed for the RandGenT
/// @param pointsPerVoxel  The number of points to scatter per voxel
/// @note returns the scattered PointDataGrid
template<
    typename GridT,
    typename RandGenT,
    typename PointDataGridT>
inline typename PointDataGridT::Ptr
pointVoxelScatter(const GridT& grid,
                  const unsigned int seed,
                  const float pointsPerVoxel);

/// @brief Non uniformly scatter points per active voxel. The pointsPerVoxel value is used
///        to weight each grids cell value to compute a fixed number of points for every
///        active voxel. If the computed result is a fractional value, each voxel calculates
///        a delta value of how likely it is to contain an extra point.
///
/// @param grid            A source grid. The resulting PointDataGrid will copy this grids
///                        transform, voxelized topology and use its values to compute a
///                        target points per voxel. The grids ValueType must be convertible
///                        to a scalar value.
/// @param seed            A seed for the RandGenT
/// @param pointsPerVoxel  The number of points to scatter per voxel
/// @note returns the scattered PointDataGrid
template<
    typename GridT,
    typename RandGenT,
    typename PointDataGridT>
inline typename PointDataGridT::Ptr
pointWeightedVoxelScatter(const GridT& grid,
                          const unsigned int seed,
                          const float pointsPerVoxel);


////////////////////////////////////////


/// Forward declaration of base class
template <typename RandGenT,
          typename ValueType,
          typename CodecT,
          typename PointDataGridT,
          typename InterruptT>
class BasePointScatter;

/// @brief The point scatters depend on the following class:
///
/// The @c InterruptType template argument below refers to any class
/// with the following interface:
/// @code
/// class Interrupter {
///   ...
/// public:
///   void start(const char* name = NULL)// called when computations begin
///   void end()                         // called when computations end
///   bool wasInterrupted(int percent=-1)// return true to break computation
///};
/// @endcode
///
/// @note If no template argument is provided for this InterruptType
/// the util::NullInterrupter is used which implies that all
/// interrupter calls are no-ops (i.e. incurs no computational overhead).

/// @brief Scatters a fixed (and integer) number of points in all
/// active voxels. The initial points per voxel value provided can be
/// a floating point value, in which case each voxel computes a delta
/// value of how likely it is to contain an extra point.
template <typename RandGenT,
          typename ValueType = float,
          typename CodecT = NullCodec,
          typename PointDataGridT = PointDataGrid,
          typename InterruptT = util::NullInterrupter>
class UniformVoxelPointScatter :
    public BasePointScatter<RandGenT, ValueType, CodecT, PointDataGridT, InterruptT>
{
public:
    using BaseT = BasePointScatter<RandGenT, ValueType, CodecT, PointDataGridT, InterruptT>;
    using BaseT::operator();

    UniformVoxelPointScatter(const RandGenT& randGen,
                             const float pointsPerVoxel,
                             const Index64 seed,
                             const float spread = 1.0f,
                             InterruptT* interrupter = nullptr)
        : BaseT(randGen, seed, spread, interrupter)
        , mPointsPerVoxel(math::Floor(pointsPerVoxel))
        , mDelta(pointsPerVoxel - float(mPointsPerVoxel))
        , mFractional(!math::isApproxZero(mDelta, 1.0e-6)) {}

    /// @brief Execute uniform voxel point scattering for a given grid of any type
    /// @note only the grids topology is used
    /// @parm grid   The source grid to base the topology generation from
    template<typename GridT>
    void operator()(const GridT& grid)
    {
        BaseT::start("Dense uniform scattering with fixed point count");
        if (!BaseT::init(grid)) return;
        if (mPointsPerVoxel < 0) return;
        if (mPointsPerVoxel == 0 && !mFractional) return;

        typename BaseT::LeafManagerT leafManager(BaseT::mPoints->tree());
        leafManager.foreach(*this);
        BaseT::end();
    }

    /// @brief Print information about the scattered points
    /// @parm name  A name to insert into the printed info
    /// @parm os    The output stream
    void print(const std::string &name, std::ostream& os = std::cout) const
    {
        if (!BaseT::mPoints) return;
        os << "Uniformly scattered " << pointCount(BaseT::mPoints->tree())
           << " points into " << BaseT::mPoints->activeVoxelCount()
           << " active voxels in \"" << name << "\" corresponding to "
           << mPointsPerVoxel << " points per voxel." << std::endl;
    }

protected:

    /// @brief Apply voxel offset values to a given leaf node
    /// @param leaf  A Point Data leaf node to apply offsets to
    /// @parm idx    The index of the leaf node
    virtual Index32 applyOffsets(typename BaseT::LeafNodeT& leaf, const size_t) const
    {
        // @note can't use iter.setValue(offset) on point grids
        Index32 offset(0);
        if (mFractional) {
            for (auto iter = leaf.beginValueAll(); iter; ++iter) {
                if (iter.isValueOn()) {
                    offset += mPointsPerVoxel;
                    if (BaseT::getRand01() < mDelta) ++offset;
                }
                leaf.setOffsetOnly(iter.pos(), offset);
            }
        }
        else {
            for (auto iter = leaf.beginValueAll(); iter; ++iter) {
                if (iter.isValueOn()) offset += mPointsPerVoxel;
                leaf.setOffsetOnly(iter.pos(), offset);
            }
        }
        return offset;
    }

private:
    const Index32 mPointsPerVoxel;
    const double mDelta;
    const bool mFractional;
};

/// @brief Scatters a varying number of points in all active voxels.
/// The value per voxel is computed by the provided points per voxel
/// multiplied by a source scalar grid. The resulting value can be
/// a floating point value, in which case each voxel computes a delta
/// value of how likely it is to contain an extra point.
template <typename RandGenT,
          typename ValueType = float,
          typename CodecT = NullCodec,
          typename PointDataGridT = PointDataGrid,
          typename InterruptT = util::NullInterrupter>
class NonUniformVoxelPointScatter :
    public BasePointScatter<RandGenT, ValueType, CodecT, PointDataGridT, InterruptT>
{
public:
    using BaseT = BasePointScatter<RandGenT, ValueType, CodecT, PointDataGridT, InterruptT>;
    using BaseT::operator();

    NonUniformVoxelPointScatter(const RandGenT& randGen,
                                const float pointsPerVoxel,
                                const Index64 seed,
                                const float spread = 1.0f,
                                InterruptT* interrupter = nullptr)
        : BaseT(randGen, seed, spread, interrupter)
        , mPointsPerVoxel(pointsPerVoxel) {}

    /// @brief Execute non uniform point scattering for a given grid of a scalar type
    /// @parm grid   The source grid to base the topology and value generation from
    template<typename GridT>
    void operator()(const GridT& grid)
    {
        BaseT::start("Non-uniform scattering with local point density");
        if (!BaseT::init(grid)) return;
        if (mPointsPerVoxel <= 0.0f) return;

        const unsigned int seed = BaseT::mSeed;
        const auto& gen = BaseT::mRand01;
        const float points = mPointsPerVoxel;
        const auto accessor = grid.getConstAccessor();

        typename BaseT::LeafManagerT leafManager(BaseT::mPoints->tree());
        leafManager.foreach([seed, accessor, points, gen]
                            (typename BaseT::LeafNodeT& leaf, const size_t idx) {
            gen.setSeed(seed + static_cast<unsigned int>(idx));
            Index32 offset(0);
            for (auto iter = leaf.beginValueAll(); iter; ++iter) {
                if (iter.isValueOn()) {
                    double fractional =
                        accessor.getValue(iter.getCoord()) * points;
                    fractional = std::max(0.0, fractional);
                    const int count = int(fractional);
                    offset += count;
                    if (gen() < (fractional - double(count))) ++offset;
                }
                // @note can't use iter.setValue(offset) on point grids
                leaf.setOffsetOnly(iter.pos(), offset);
            }
        });

        leafManager.foreach(*this);
        BaseT::end();
    }

    /// @brief Print information about the scattered points
    /// @parm name  A name to insert into the printed info
    /// @parm os    The output stream
    void print(const std::string &name, std::ostream& os = std::cout) const
    {
        if (!BaseT::mPoints) return;
        os << "Non-uniformly scattered " << pointCount(BaseT::mPoints->tree())
           << " points into " << BaseT::mPoints->activeVoxelCount()
           << " active voxels in \"" << name << "\"." << std::endl;
    }

protected:

    /// @brief A nullop which simply returns the amount of desired points.
    /// The value offsets are generated in the above operator
    /// @param leaf  A Point Data leaf node
    /// @parm idx    The index of the leaf node
    virtual Index32 applyOffsets(typename BaseT::LeafNodeT& leaf, const size_t) const
    {
        // Offsets are applied in first foreach in the above operator
        return leaf.getLastValue();
    }

private:
    const float mPointsPerVoxel;
};

/// @brief Uniformly scatters a total integer amount of points in the
/// active voxels. Each voxel receives a base amount of points, calculated
/// from the total integer point count provided and the number of active
/// voxels. The remaining values are randomly scattered.
template <typename RandGenT,
          typename ValueType = float,
          typename CodecT = NullCodec,
          typename PointDataGridT = PointDataGrid,
          typename InterruptT = util::NullInterrupter>
class TotalPointScatter :
    public BasePointScatter<RandGenT, ValueType, CodecT, PointDataGridT, InterruptT>
{
public:
    using BaseT = BasePointScatter<RandGenT, ValueType, CodecT, PointDataGridT, InterruptT>;
    using BaseT::operator();

    TotalPointScatter(const RandGenT& randGen,
                      const Index64 totalPointCount,
                      const Index64 seed,
                      const float spread = 1.0f,
                      InterruptT* interrupter = nullptr)
        : BaseT(randGen, seed, spread, interrupter)
        , mPointCount(totalPointCount) {}

    /// @brief Execute total point scattering for a given grid of any type
    /// @note only the grids topology is used
    /// @parm grid   The source grid to base the topology generation from
    template<typename GridT>
    void operator()(const GridT& grid)
    {
        BaseT::start("Uniform scattering with fixed point count");
        if (!BaseT::init(grid)) return;
        if (mPointCount == 0) return;

        typename BaseT::LeafManagerT leafManager(BaseT::mPoints->tree());
        const Index64 voxelCount = leafManager.activeLeafVoxelCount();
        assert(voxelCount != 0);

        const double pointsPerVolume = double(mPointCount) / double(voxelCount);
        mPointsPerVoxel = static_cast<Index32>(math::RoundDown(pointsPerVolume));
        Index32 remainder = mPointCount - (mPointsPerVoxel * voxelCount);

        if (remainder == 0) {
            UniformVoxelPointScatter<RandGenT, ValueType, CodecT, PointDataGridT, InterruptT>
                op(BaseT::mRand01.engine(), float(mPointsPerVoxel),
                 BaseT::mSeed, BaseT::mSpread, BaseT::mInterrupter);
            op(grid);
            BaseT::mPoints = op.points();
            return;
        }

        const auto leafRange = leafManager.leafRange();

        Index64 offset = 0;
        std::vector<Index64> voxelOffsets;
        voxelOffsets.reserve(leafManager.leafCount());
        for (auto leaf = leafRange.begin(); leaf; ++leaf) {
            offset += leaf->onVoxelCount();
            voxelOffsets.push_back(offset);
        }

        std::vector<Index64> values(remainder);
        const unsigned int seed = BaseT::mSeed;
        using RangeT = tbb::blocked_range<size_t>;
        tbb::parallel_for(RangeT(0, remainder), [voxelCount, &values, seed]
                          (const RangeT& range) {
            size_t n = range.begin();
            math::RandInt<Index64, RandGenT>
                gen(static_cast<unsigned int>(n) + seed, 0, voxelCount-1);
            for (size_t N = range.end(); n < N; ++n) {
                values[n] = gen();
            }
        });
        tbb::parallel_sort(values.begin(), values.end());

        mVoxelOffsets = &voxelOffsets;
        mValues = &values;
        leafManager.foreach(*this);
        BaseT::end();
    }

    /// @brief Print information about the scattered points
    /// @parm name  A name to insert into the printed info
    /// @parm os    The output stream
    void print(const std::string &name, std::ostream& os = std::cout) const
    {
        if (!BaseT::mPoints) return;
        const Index64 points = pointCount(BaseT::mPoints->tree());
        const Index64 voxels = BaseT::mPoints->activeVoxelCount();

        os << "Uniformly scattered " << points << " points into " << voxels
           << " active voxels in \"" << name << "\" corresponding to "
           << (double(points) / double(voxels)) << " points per voxel." << std::endl;
    }

protected:

    /// @brief Apply voxel offset values to a given leaf node
    /// @param leaf  A Point Data leaf node to apply offsets to
    /// @parm idx    The index of the leaf node
    virtual Index32 applyOffsets(typename BaseT::LeafNodeT& leaf, const size_t idx) const
    {
        const Index64 lowerOffset = idx == 0 ? 0 : (*mVoxelOffsets)[idx-1];
        const Index64 vOffset = (*mVoxelOffsets)[idx];
        assert(vOffset > lowerOffset);

        const auto valuesBegin = mValues->begin();
        auto lower = std::lower_bound(valuesBegin, mValues->end(), vOffset);
        assert(lower != valuesBegin);

        auto* const data = leaf.buffer().data();
        while (lower-- != valuesBegin) {
            const Index64 vId = *lower;
            if (vId < lowerOffset) break;
            auto& offset = data[vId - lowerOffset];
            offset = offset + 1; // no += operator support
        }

        Index32 offset(0);
        for (auto iter = leaf.beginValueAll(); iter; ++iter) {
            if (iter.isValueOn()) offset += mPointsPerVoxel + data[iter.pos()];
            // @note can't use iter.setValue(offset) on point grids
            leaf.setOffsetOnly(iter.pos(), offset);
        }
        return offset;
    }

private:
    const Index64 mPointCount;
    Index32 mPointsPerVoxel;
    std::vector<Index64>* mVoxelOffsets;
    std::vector<Index64>* mValues;
};


template <typename RandGenT,
          typename ValueType,
          typename CodecT,
          typename PointDataGridT,
          typename InterruptT>
class BasePointScatter
{
public:
    using RandomGenerator = math::Rand01<ValueType, RandGenT>;
    using TreeT = typename PointDataGridT::TreeType;
    using LeafNodeT = typename TreeT::LeafNodeType;
    using LeafManagerT = tree::LeafManager<TreeT>;

    /// @brief Retrieve the scattered points
    /// @note will contain a nullptr if the derived operator
    /// has not been called
    inline typename PointDataGridT::Ptr points()
    {
        return mPoints;
    }

    /// @brief initialise the RNG seed
    /// @parm seed   The seed for the random number generator
    inline void setSeed(const unsigned int seed)
    {
        mSeed = seed;
    }

    /// @brief set the spread value
    /// @parm spread   The spread of points from the voxel center from 0 to 1
    inline void setSpread(const float spread)
    {
        assert(spread >= 0.0f && spread <= 1.0f);
        mSpread = spread;
    }

    /// @brief The main operator for building leaf node positions. This is
    /// called by the derived scatting classes through a parallel executer
    /// @param leaf  A Point Data leaf node to scatter
    /// @parm idx    The index of the leaf node
    void operator()(LeafNodeT& leaf, const size_t idx) const
    {
        assert(!leaf.empty());
        if (util::wasInterrupted(mInterrupter)) {
            tbb::task::self().cancel_group_execution();
            return;
        }
        mRand01.setSeed(static_cast<unsigned int>(idx) + mSeed);
        const Index32 count = this->applyOffsets(leaf, idx);
        this->initializePositions(leaf, count);
    }

protected:

    /// @brief The virtual method to compute point offsets for each leaf which
    /// derived classes must provide. This is called for each leaf node in
    /// generated the PointDataTree topology
    /// @param leaf  A Point Data leaf node to apply offsets to
    /// @parm idx    The index of the leaf node
    virtual Index32 applyOffsets(LeafNodeT& leaf, const size_t idx) const = 0;

    /// @brief initialize the interrupter message
    /// @parm name   A message to provide to the interrupter
    inline void start(const char* name)
    {
        if (mInterrupter) mInterrupter->start(name);
    }

    /// @brief end the interrupter message
    inline void end()
    {
        if (mInterrupter) mInterrupter->end();
    }

    /// @brief initialise the topology of a PointDataGrid and ensure
    /// everything is voxelized
    /// @parm grid   The source grid to base the topology generation from
    template<typename GridT>
    inline bool init(const GridT& grid)
    {
        mPoints.reset(new PointDataGridT);
        mPoints->setTransform(grid.transform().copy());
        mPoints->topologyUnion(grid);
        if (mPoints->tree().hasActiveTiles()) {
            mPoints->tree().voxelizeActiveTiles();
        }

        return static_cast<bool>(mPoints->constTree().cbeginLeaf());
    }

    /// @brief Return a random floating point number between zero and one
    inline double getRand01() const
    {
        return mRand01();
    }

    /// @brief Return a random floating point number between -0.5 -+ mSpread/2
    inline ValueType getRand() const
    {
        return mSpread * (mRand01() - ValueType(0.5));
    }

private:

    using PositionValueT = math::Vec3<ValueType>;
    using PositionArrayT = TypedAttributeArray<PositionValueT, CodecT>;
    using PositionWriteHandle = AttributeWriteHandle<PositionValueT, CodecT>;

    /// @brief Generate a new points voxel position
    /// @parm pos   The position to set
    inline void getPoint(PositionValueT& pos) const
    {
        pos[0] = getRand();
        pos[1] = getRand();
        pos[2] = getRand();
    }

    /// @brief Populate the voxel position array for a leaf node
    /// @parm leaf   The leaf node to build positions for
    /// @parm count  The total number of points to allocate
    inline void initializePositions(LeafNodeT& leaf, const Index64 count) const
    {
        leaf.initializeAttributes(mAttributeDescriptor, count);
        PositionWriteHandle pHandle(leaf.attributeArray(0));
        PositionValueT P;
        for (Index64 index = 0; index < count; ++index) {
            getPoint(P);
            pHandle.set(index, P);
        }
    }

protected:

    /// @brief Base constructor
    /// @param randGen      The random number generator to use
    /// @param seed         The seed for the random number generator
    /// @param spread       The spread of voxel positions from the center of a voxel
    /// @param interrupter  An optional interrupter
    BasePointScatter(const RandGenT& randGen,
                     const unsigned int seed,
                     const float spread,
                     InterruptT* interrupter)
        : mPoints(nullptr)
        , mRand01(randGen)
        , mSpread(spread)
        , mSeed(seed)
        , mInterrupter(interrupter)
        , mAttributeDescriptor(AttributeSet::Descriptor::create(PositionArrayT::attributeType()))
        { }

    typename PointDataGridT::Ptr mPoints;
    mutable RandomGenerator mRand01;
    ValueType mSpread;
    unsigned int mSeed;
    InterruptT* mInterrupter;

private:
    AttributeSet::Descriptor::Ptr mAttributeDescriptor;
};


template<
    typename GridT,
    typename RandGenT,
    typename PointDataGridT>
inline typename PointDataGridT::Ptr
pointScatter(const GridT& grid,
             const unsigned int seed,
             const Index64 count)
{
    using ScatterT = TotalPointScatter<RandGenT, float,
        NullCodec, PointDataGridT>;

    RandGenT gen(seed);
    ScatterT scatter(gen, count, seed);
    scatter(grid);
    return scatter.points();
}

template<
    typename GridT,
    typename RandGenT,
    typename PointDataGridT>
inline typename PointDataGridT::Ptr
pointVoxelScatter(const GridT& grid,
                  const unsigned int seed,
                  const float pointsPerVoxel)
{
    using ScatterT = UniformVoxelPointScatter<RandGenT, float,
        NullCodec, PointDataGridT>;

    RandGenT gen(seed);
    ScatterT scatter(gen, pointsPerVoxel, seed);
    scatter(grid);
    return scatter.points();
}

template<
    typename GridT,
    typename RandGenT,
    typename PointDataGridT>
inline typename PointDataGridT::Ptr
pointWeightedVoxelScatter(const GridT& grid,
                          const unsigned int seed,
                          const float pointsPerVoxel)
{
    using ScatterT = NonUniformVoxelPointScatter<RandGenT, float,
        NullCodec, PointDataGridT>;

    RandGenT gen(seed);
    ScatterT scatter(gen, pointsPerVoxel, seed);
    scatter(grid);
    return scatter.points();
}


} // namespace points
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb


#endif // OPENVDB_POINTS_POINT_SCATTER_HAS_BEEN_INCLUDED

// Copyright (c) 2012-2017 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
